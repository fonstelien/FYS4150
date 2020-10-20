#include "catch.hpp"
#include "utils.hpp"

using namespace arma;
using namespace std;

#define TOL 1.E-3

/* From utils.cpp */


void get_pos_and_vel(int i, mat &flight_log, vec &pos, vec &vel) {
  pos(0) = flight_log(i,0);
  pos(1) = flight_log(i,1);

  vel(0) = flight_log(i,2);
  vel(1) = flight_log(i,3);
}


TEST_CASE("Verlet energy and position tests", "[verlet-tests]") {
  int n = 10000;
  double h = 1./n;
  double ep0, ep, ek0, ek;
  mat flight_log;
  vec pos0(2), pos(2), vel0(2), vel(2);
  
  flight_log = earth_circular_verlet(h, n+1);
  get_pos_and_vel(0, flight_log, pos0, vel0);
  get_pos_and_vel(n, flight_log, pos, vel);

  // Test preservation of potential and kinetic energy
  ep0 = norm(pos0);
  ep = norm(pos);
  ek0 = norm(vel0)*norm(vel0);
  ek = norm(vel)*norm(vel);
  REQUIRE(ep == Approx(ep0).epsilon(1E-6));
  REQUIRE(ek == Approx(ek0).epsilon(1E-6));

  // Test arrival at initial position
  REQUIRE(pos(0) == Approx(pos0(0)).margin(1E-6));
  REQUIRE(pos(1) == Approx(pos0(1)).margin(1E-6));
}


TEST_CASE("Euler energy and position tests", "[euler-tests]") {
  int n = 400000;
  double h = 1./n;
  double ep0, ep, ek0, ek;
  mat flight_log;
  vec pos0(2), pos(2), vel0(2), vel(2);
  
  flight_log = earth_circular_fwd_euler(h, n+1);
  get_pos_and_vel(0, flight_log, pos0, vel0);
  get_pos_and_vel(n, flight_log, pos, vel);

  // Test preservation of potential and kinetic energy
  ep0 = norm(pos0);
  ep = norm(pos);
  ek0 = norm(vel0)*norm(vel0);
  ek = norm(vel)*norm(vel);
  REQUIRE(ep == Approx(ep0).epsilon(1.0E-3));
  REQUIRE(ek == Approx(ek0).epsilon(1.0E-3));

  // Test arrival at initial position
  REQUIRE(pos(0) == Approx(pos0(0)).margin(1.0E-3));
  REQUIRE(pos(1) == Approx(pos0(1)).margin(1.0E-3));
}


TEST_CASE("System potential, kinetic energy; angular momentum", "[sys-pot-kin-mom]") {
  int n = 100000;
  double h = 1./n;
  vec pos(3), vel(3);
  double ep0, ep, ek0, ek, am0, am;

  pos = {0.,0.,0.};
  vel = {0.,0.,0.};
  Planet sun = Planet(1., pos, vel);
  sun.fixed = true;
  
  pos = {1.,0.,0.};
  vel = {0.,5.,0.};  // elliptic orbit
  Planet earth = Planet(3.0E-6, pos, vel);

  Solver solver;
  solver.add(&sun);
  solver.add(&earth);

  ep0 = solver.potential_energy(&earth, &sun);
  ek0 = earth.kinetic_energy();
  am0 = earth.angular_momentum();

  solver.solve(n/3, h);
  ep = solver.potential_energy(&earth, &sun);
  ek = earth.kinetic_energy();
  am = earth.angular_momentum();
  REQUIRE(ep+ek == Approx(ep0+ek0).epsilon(1.0E-3));
  REQUIRE(am == Approx(am0).epsilon(1.0E-6));

  solver.solve(n/3, h);
  ep = solver.potential_energy(&earth, &sun);
  ek = earth.kinetic_energy();
  am = earth.angular_momentum();
  REQUIRE(ep+ek == Approx(ep0+ek0).epsilon(1.0E-3));
  REQUIRE(am == Approx(am0).epsilon(1.0E-6));

  solver.solve(n/3, h);
  ep = solver.potential_energy(&earth, &sun);
  ek = earth.kinetic_energy();
  am = earth.angular_momentum();
  REQUIRE(ep+ek == Approx(ep0+ek0).epsilon(1.0E-3));
  REQUIRE(am == Approx(am0).epsilon(1.0E-6));
}


TEST_CASE("Forms of the force", "[forms-of-the-force]") {
  int n = 100000;
  double h = 1./n;
  vec pos(3), vel(3);
  double ep0, ep, ek0, ek, am0, am;

  pos = {0.,0.,0.};
  vel = {0.,0.,0.};
  Planet sun = Planet(1., pos, vel);
  sun.fixed = true;
  
  pos = {1.,0.,0.};
  vel = {0.,5.,0.};  // elliptic orbit
  Planet earth = Planet(3.0E-6, pos, vel);

  Solver solver;
  solver.add(&sun);
  solver.add(&earth);

  n = (int) n * .26;
  solver.beta = 3.0;
  
  ep0 = solver.potential_energy(&earth, &sun);
  ek0 = earth.kinetic_energy();
  am0 = earth.angular_momentum();

  solver.solve(n/3, h);
  ep = solver.potential_energy(&earth, &sun);
  ek = earth.kinetic_energy();
  am = earth.angular_momentum();
  REQUIRE(ep+ek == Approx(ep0+ek0).epsilon(1.0E-3));
  REQUIRE(am == Approx(am0).epsilon(1.0E-6));

  solver.solve(n/3, h);
  ep = solver.potential_energy(&earth, &sun);
  ek = earth.kinetic_energy();
  am = earth.angular_momentum();
  REQUIRE(ep+ek == Approx(ep0+ek0).epsilon(1.0E-3));
  REQUIRE(am == Approx(am0).epsilon(1.0E-6));

  solver.solve(n/3, h);
  ep = solver.potential_energy(&earth, &sun);
  ek = earth.kinetic_energy();
  am = earth.angular_momentum();
  REQUIRE(ep+ek == Approx(ep0+ek0).epsilon(1.0E-3));
  REQUIRE(am == Approx(am0).epsilon(1.0E-6));
}
