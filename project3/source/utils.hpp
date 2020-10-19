#include <iostream>
#include <armadillo>
#include <cmath>
#include <vector>

#define DEBUG(msg) cout << "DEBUG " << msg << endl;

// four-pi-squared
#define FPS 39.47841760435743

// Crashed planet threshold
#define CRASH_LIM 1.0E-6

// Speed of light in AU per year
#define LIGHTPSPEED (2.99792458E+8 / 1.495978707E+11 * 60*60*24*365)

using namespace arma;
using namespace std;

/* From utils.cpp */

/* From non_OO_earth_sun.cpp */
mat earth_circular_fwd_euler(double h, int n);
mat earth_circular_verlet(double h, int n);
mat time_algorithms();

/* From perihelion_of_mercury.cpp */
double perihelion_of_mercury(double h, long long int n, bool general_relativistic);

/* From planet.cpp */
class Planet {
private:
  vec acc0, acc;  // acceleration x,y,z directions

public:
  double m;  // mass relative to the sun
  vec pos, vel;  // position, velocity x,y,z directions
  bool fixed = false;  // planet has fixed position
  bool crashed = false;  // planet has crashed
  
  // Constructor
  Planet(double mass, vec init_pos, vec init_vel);

  // Velocity Verlet functions
  // Call update_acc() once to initiate, then in this order:
  void update_pos(double dt);  // update position
  void update_acc(Planet *other, vec r, double r3);  // update acceleration towards other
  void update_vel(double dt);  // update velocity
  double kinetic_energy();
  double angular_momentum();
};


/* From solver.cpp */
class Solver {
public:
  vector<Planet *> planets;  // collection of planets
  mat flight_log;  // log of trajectories and velocities of all planets in colleciton
  double beta = 2.;  // experimental exponent in Newton's law of gravitation
  
  void add(Planet *planet);  // add planet to collection
  void build(mat system);  // builds system from system
  void nuke();  // frees Planets pointed to by planets and clear()s the vector
  void solve(int steps, double dt);  // run simulation with num. steps, time step dt
  string csv_header();  // returns csv format header for the flight_log
  double potential_energy(Planet *p1, Planet *p2);  // potential energy between p1 and p2
};
