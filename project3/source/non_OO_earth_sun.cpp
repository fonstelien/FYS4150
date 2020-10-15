#include "utils.hpp"

using namespace arma;
using namespace std;

// Earth's initial position
#define X0 1.
#define Y0 0.

void euler_fwd_iteration(double h, vec &pos, vec &vel) {
  vec pos0 = pos;
  vec vel0 = vel;
  
  pos = pos0 + h*vel0;
  vel = vel0 - h*FPS/norm(pos0)*pos0;
}

void velocity_verlet_iteration(double h, vec &pos, vec &vel) {
  vec pos0 = pos;
  vec vel0 = vel;
  vec acc0(3), acc1(3);

  acc0 = -FPS/norm(pos0)*pos0;
  pos = pos0 + h*vel0 + h*h/2*acc0;
  acc1 = -FPS/norm(pos)*pos;  
  vel = vel0 + h/2*(acc1 + acc0);
}
