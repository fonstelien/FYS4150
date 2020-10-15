#include "utils.hpp"

using namespace arma;
using namespace std;


/* Set Earth's initial position and velocity */
void initialize_earth_circular_trajectory(vec &pos, vec &vel) {
  pos(0) = 1.;
  pos(1) = 0.;
  pos(2) = 0.;
  vel(0) = 0.;
  vel(1) = 2*M_PI;  // 1 rot per annum
  vel(2) = 0.;
}


/* Runs forward euler with earth in circular trajectory */
void earth_circular_fwd_euler(double h, mat &results) {
  int n = results.n_rows;
  vec pos0(3), pos(3), vel(3);
  double tmp;
  
  initialize_earth_circular_trajectory(pos, vel);
  results(0,0) = pos(0);
  results(0,1) = pos(1);
  results(0,2) = pos(2);
  results(0,3) = vel(0);
  results(0,4) = vel(1);
  results(0,5) = vel(2);

  for (int i = 1; i < n; i++) {
    pos0 = pos;
    pos += h*vel;
    tmp = norm(pos0);
    vel -= h*FPS/(tmp*tmp*tmp)*pos0;
    
    results(i,0) = pos(0);
    results(i,1) = pos(1);
    results(i,2) = pos(2);
    results(i,3) = vel(0);
    results(i,4) = vel(1);
    results(i,5) = vel(2);
  }
}


/* Runs verlet with earth in circular trajectory */
void earth_circular_verlet(double h, mat &results) {
  int n = results.n_rows;
  vec pos(3), vel(3), acc0(3), acc1(1);
  double tmp;
  
  initialize_earth_circular_trajectory(pos, vel);
  results(0,0) = pos(0);
  results(0,1) = pos(1);
  results(0,2) = pos(2);
  results(0,3) = vel(0);
  results(0,4) = vel(1);
  results(0,5) = vel(2);
  
  for (int i = 0; i < n; i++) {
    tmp = norm(pos);
    acc0 = -FPS/(tmp*tmp*tmp)*pos;
    pos += h*vel + h*h/2*acc0;
    tmp = norm(pos);
    acc1 = -FPS/(tmp*tmp*tmp)*pos;  
    vel += h/2*(acc1 + acc0);
    
    results(i,0) = pos(0);
    results(i,1) = pos(1);
    results(i,2) = pos(2);
    results(i,3) = vel(0);
    results(i,4) = vel(1);
    results(i,5) = vel(2);
  }
}



/* One iteration with the forward Euler algorithm */
void euler_fwd_iteration(double h, vec &pos, vec &vel) {
  vec pos0 = pos;
  vec vel0 = vel;
  
  pos = pos0 + h*vel0;
  vel = vel0 - h*FPS/norm(pos0)*pos0;
}


/* One iteration with the forward Velocity Verlet algorithm */
void velocity_verlet_iteration(double h, vec &pos, vec &vel) {
  vec pos0 = pos;
  vec vel0 = vel;
  vec acc0(3), acc1(3);

  acc0 = -FPS/norm(pos0)*pos0;
  pos = pos0 + h*vel0 + h*h/2*acc0;
  acc1 = -FPS/norm(pos)*pos;  
  vel = vel0 + h/2*(acc1 + acc0);
}
