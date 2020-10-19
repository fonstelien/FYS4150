#include <time.h>
#include "utils.hpp"

using namespace arma;
using namespace std;


/* Set Earth's initial position and velocity */
void initialize_earth_circular_trajectory(vec &pos, vec &vel) {
  pos(0) = 1.;  // distance 1 AU from the sun
  pos(1) = 0.;
  pos(2) = 0.;
  vel(0) = 0.;
  vel(1) = 2*M_PI;  // 1 rot per annum
  vel(2) = 0.;
}


/* Runs forward euler with earth in circular orbit and logs */
/* trajectory and velocity in flight_log. */
mat earth_circular_fwd_euler(double h, int n) {
  mat flight_log(n,6);
  vec pos0(3), pos(3), vel(3);
  double tmp;
  
  initialize_earth_circular_trajectory(pos, vel);
  flight_log(0,0) = pos(0);
  flight_log(0,1) = pos(1);
  flight_log(0,2) = pos(2);
  flight_log(0,3) = vel(0);
  flight_log(0,4) = vel(1);
  flight_log(0,5) = vel(2);

  for (int i = 1; i < n; i++) {
    pos0 = pos;
    pos += h*vel;
    tmp = norm(pos0);
    vel -= h*FPS/(tmp*tmp*tmp)*pos0;
    
    flight_log(i,0) = pos(0);
    flight_log(i,1) = pos(1);
    flight_log(i,2) = pos(2);
    flight_log(i,3) = vel(0);
    flight_log(i,4) = vel(1);
    flight_log(i,5) = vel(2);
  }

  return flight_log;  
}


/* Runs Velocity Verlet with earth in circular orbit and logs */
/* trajectory and velocity in flight_log. */
mat earth_circular_verlet(double h, int n) {
  mat flight_log(n,6);
  vec pos(3), vel(3), acc0(3), acc(3);
  double tmp;
  
  initialize_earth_circular_trajectory(pos, vel);
  flight_log(0,0) = pos(0);
  flight_log(0,1) = pos(1);
  flight_log(0,2) = pos(2);
  flight_log(0,3) = vel(0);
  flight_log(0,4) = vel(1);
  flight_log(0,5) = vel(2);

  // Initializing
  tmp = norm(pos);
  acc = -FPS/(tmp*tmp*tmp)*pos;

  for (int i = 1; i < n; i++) {
    pos += h*vel + h*h/2*acc;
    acc0 = acc;
    tmp = norm(pos);
    acc = -FPS/(tmp*tmp*tmp)*pos;  
    vel += h/2*(acc + acc0);
    
    flight_log(i,0) = pos(0);
    flight_log(i,1) = pos(1);
    flight_log(i,2) = pos(2);
    flight_log(i,3) = vel(0);
    flight_log(i,4) = vel(1);
    flight_log(i,5) = vel(2);
  }

  return flight_log;
}


mat time_algorithms() {
  clock_t start, stop;
  double lapsed_time;
  int T = 9;
  vec N = logspace(1,T,T);
  mat time_log(T,3);
  int runs = 10;
  int n;
  double h;
  vec pos0(3), pos(3), vel(3), acc0(3), acc(3);
  double tmp;

  // Time Euler runs
  for (int i = 0; i < T; i++) {
    n = N(i);
    h = 1./n;
    start = clock();
    for (int r = 0; r < runs; r++) {
      initialize_earth_circular_trajectory(pos, vel);
      for (int j = 0; j < n; j++) {
	pos0 = pos;
	pos += h*vel;
	tmp = norm(pos0);
	vel -= h*FPS/(tmp*tmp*tmp)*pos0;
      }
    }
    stop = clock();
    lapsed_time = (double) (stop-start) / CLOCKS_PER_SEC;
    time_log(i,0) = (double) n;
    time_log(i,1) = lapsed_time / runs;
  }

  // Time Verlet runs
  for (int i = 0; i < T; i++) {
    n = N(i);
    h = 1./n;
    start = clock();
    for (int r = 0; r < runs; r++) {
      initialize_earth_circular_trajectory(pos, vel);
      tmp = norm(pos);
      acc = -FPS/(tmp*tmp*tmp)*pos;
      for (int j = 0; j < n; j++) {
	pos += h*vel + h*h/2*acc;
	acc0 = acc;
	tmp = norm(pos);
	acc = -FPS/(tmp*tmp*tmp)*pos;  
	vel += h/2*(acc + acc0);

      }
      stop = clock();
      lapsed_time = (double) (stop-start) / CLOCKS_PER_SEC;
      time_log(i,2) = lapsed_time / runs;
    }
  }
  
  return time_log;
}
