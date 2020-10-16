#include <iostream>
#include <armadillo>
#include <cmath>
#include <vector>

#define DEBUG(msg) cout << "DEBUG " << msg << endl;

// four-pi-squared
#define FPS 39.47841760435743

// Mass of the sun
#define M_SUN 2.0E30

using namespace arma;
using namespace std;

/* From utils.cpp */

/* From non_OO_earth_sun.cpp */
void earth_circular_fwd_euler(double h, mat &results);
void earth_circular_verlet(double h, mat &results);


/* From planet.cpp */
class Planet {
private:
  vec acc0, acc;  // acceleration x,y,z directions

public:
  double m;  // mass relative to the sun
  vec pos, vel;  // position, velocity x,y,z directions
  bool fixed = false;  // fixed position
  
  // Constructor
  Planet(double mass, vec init_pos, vec init_vel);

  // Velocity Verlet functions
  // Call update_acc() once to initiate, then in this order:
  void update_pos(double dt);  // update position
  void update_acc(Planet *other, vec r, double r3);  // update acceleration towards other
  void update_vel(double dt);  // update velocity
};


/* From solver.cpp */
class Solver {
public:
  vector<Planet> planets;  // collection of planets
  mat flight_log;  // log of trajectories and velocities of all planets in colleciton

  void add(Planet planet);  // add planet to collection
  void build(mat system);  // builds system from system
  void solve(int steps, double dt);  // run simulation with num. steps, time step dt
  string csv_header();  // returns csv format header for the flight_log
};
