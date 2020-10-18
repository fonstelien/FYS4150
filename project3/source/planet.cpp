#include "utils.hpp"

using namespace arma;
using namespace std;


/* Constructor */
/* Mass mass relative to sun; initial position and velocity vecs */
Planet::Planet(double mass, vec init_pos, vec init_vel) {
  m = mass;
  pos = init_pos;
  vel = init_vel;

  acc0 = {0.,0.,0.};
  acc = {0.,0.,0.};
}


/* Updates planet's position with time step dt */
void Planet::update_pos(double dt) {
  if (fixed || crashed)
    return;

  pos += dt*vel + dt*dt/2*acc;
  acc0 = acc;
  acc = {0.,0.,0.};
}


/* Updates planet's acceleration towards other Planet in direction r, r3 is dist^3 */
void Planet::update_acc(Planet *other, vec r, double r3) {
  if (fixed || crashed)
    return;

  if (r3 <= CRASH_LIM) {
    crashed = true;
    return;
  }

  acc += FPS*other->m/r3*r;
}


/* Updates planet's velocity with time step dt */
void Planet::update_vel(double dt) {
  if (fixed || crashed)
    return;

  vel += dt/2*(acc + acc0);
}


/* Kinetic energy of this planet in sun-masses*AU^2/yr.^2 */
double Planet::kinetic_energy() {
  double v = norm(vel);
  return .5*m*v*v;
}


/* Magnitude of angular momentum of this planet about x=0,y=0,z=0 in sun-masses*AU^2/yr. */
double Planet::angular_momentum() {
  vec L = cross(pos, m*vel);
  return norm(L);
}
