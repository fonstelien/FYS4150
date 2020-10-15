#include "utils.hpp"

using namespace arma;
using namespace std;

Planet::Planet(double mass, vec init_pos, vec init_vel) {
  m = mass;
  pos = init_pos;
  vel = init_vel;

  acc0 = {0.,0.,0.};
  acc1 = {0.,0.,0.};
}

void Planet::update_acc0(Planet *other, vec r, double r3) {
  acc0 += FPS*other->m/r3*r;
}

void Planet::update_acc1(Planet *other, vec r, double r3) {
  acc1 += FPS*other->m/r3*r;
}

void Planet::update_pos(double dt) {
  pos += dt*vel + dt*dt/2*acc0;
}

void Planet::update_vel(double dt) {
  vel += dt/2*(acc1 + acc0);
  acc0 = {0.,0.,0.};
  acc1 = {0.,0.,0.};
}
