#include <iostream>
#include <armadillo>
#include <cmath>

#define DEBUG(msg) cout << "DEBUG " << msg << endl;

// four-pi-squared
#define FPS 39.47841760435743


using namespace arma;
using namespace std;

/* From utils.cpp */



/* From non_OO_earth_sun.cpp */
void euler_fwd_iteration(double h, vec &pos, vec &vel);
void velocity_verlet_iteration(double h, vec &pos, vec &vel);
