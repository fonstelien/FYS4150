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
  vec acc0, acc1;  // acceleration x,y,z directions

public:
  double m;  // mass relative to the sun
  vec pos, vel;  // position, velocity x,y,z directions

  // Constructor
  Planet(double mass, vec init_pos, vec init_vel);

  void update_acc0(Planet *other, vec r, double r3);
  void update_acc1(Planet *other, vec r, double r3);
  void update_pos(double dt);
  void update_vel(double dt);
};

/* From solver.cpp */
class Solver {
public:
  vector<Planet> planets;
  mat flight_log;

  void add(Planet planet);
  void solve(int steps, double dt);
  void to_csv();
};
