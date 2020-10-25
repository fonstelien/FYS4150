#include <time.h>
#include "utils.hpp"

using namespace arma;
using namespace std;


/* Set Earth's initial position and velocity */
void initialize_earth_circular_trajectory(vec &pos, vec &vel) {
  pos(0) = 1.;  // distance 1 AU from the sun
  pos(1) = 0.;
  vel(0) = 0.;
  vel(1) = 2*M_PI;  // one circular rotation per annum
}


/* Runs forward euler with earth in circular orbit and logs */
/* trajectory and velocity in flight_log. Calculations in xy plane only. */
mat earth_circular_fwd_euler(double h, int n) {
  mat flight_log(n+1,4);
  vec pos0(2), pos(2), vel(2);
  double r;
  
  initialize_earth_circular_trajectory(pos, vel);
  flight_log(0,0) = pos(0);
  flight_log(0,1) = pos(1);
  flight_log(0,2) = vel(0);
  flight_log(0,3) = vel(1);

  for (int i = 1; i < n+1; i++) {
    pos0 = pos;
    pos += h*vel;
    r = norm(pos0);
    vel -= h*FPS/(r*r*r)*pos0;
    
    flight_log(i,0) = pos(0);
    flight_log(i,1) = pos(1);
    flight_log(i,2) = vel(0);
    flight_log(i,3) = vel(1);
  }

  return flight_log;  
}


/* Runs Velocity Verlet with earth in circular orbit and logs */
/* trajectory and velocity in flight_log. Calculations in xy plane only. */
mat earth_circular_verlet(double h, int n) {
  mat flight_log(n+1,4);
  vec pos(2), vel(2), acc0(2), acc(2);
  double r;
  
  initialize_earth_circular_trajectory(pos, vel);
  flight_log(0,0) = pos(0);
  flight_log(0,1) = pos(1);
  flight_log(0,2) = vel(0);
  flight_log(0,3) = vel(1);

  // Initializing
  r = norm(pos);
  acc = -FPS/(r*r*r)*pos;

  for (int i = 1; i < n+1; i++) {
    pos += h*vel + h*h/2*acc;
    acc0 = acc;
    r = norm(pos);
    acc = -FPS/(r*r*r)*pos;  
    vel += h/2*(acc + acc0);
    
    flight_log(i,0) = pos(0);
    flight_log(i,1) = pos(1);
    flight_log(i,2) = vel(0);
    flight_log(i,3) = vel(1);
  }

  return flight_log;
}


/* Runs both Euler and Verlet algorithms and times them. */
/* Reported times are averaged over 10 runs. */
/* NOTE: Even if the orbit is circular in the xy plane, I still want the */
/* time for three-dimensional calculations for the report. */
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
  double r;

  // Time Euler runs
  for (int i = 0; i < T; i++) {
    n = N(i);
    h = 1./n;
    start = clock();
    for (int run = 0; run < runs; run++) {
      initialize_earth_circular_trajectory(pos, vel);
      pos(2) = 0.;
      vel(2) = 0.;
      for (int j = 0; j < n; j++) {
	pos0 = pos;
	pos += h*vel;
	r = norm(pos0);
	vel -= h*FPS/(r*r*r)*pos0;
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
    for (int run = 0; run < runs; run++) {
      initialize_earth_circular_trajectory(pos, vel);
      pos(2) = 0.;
      vel(2) = 0.;
      r = norm(pos);
      acc = -FPS/(r*r*r)*pos;
      for (int j = 0; j < n; j++) {
	pos += h*vel + h*h/2*acc;
	acc0 = acc;
	r = norm(pos);
	acc = -FPS/(r*r*r)*pos;  
	vel += h/2*(acc + acc0);
      }
      stop = clock();
      lapsed_time = (double) (stop-start) / CLOCKS_PER_SEC;
      time_log(i,2) = lapsed_time / runs;
    }
  }
  
  return time_log;
}
