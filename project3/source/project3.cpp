#include <fstream>
#include <unistd.h>
#include <time.h>
#include "utils.hpp"

/* Defaults */

/* Modes */
#define EULER 1
#define VERLET 2
#define TIME 3
#define CREATOR 4
#define PERIHELION 5
#define SYSTEM 9

/* Algorithms */

using namespace arma;
using namespace std;

/* Prints program usage */
void print_usage() {
  cout << "usage: project3 [-h | -f | -E | -V | -T | -C | -P] -n -y -b -v -j -s -r" << endl;
  cout << "  -h   Print this text" << endl;
  cout << "  -f   Run system from csv file input. Format: mass,x,y,z,vx,vy,vz" << endl;
  cout << "       mass is planet's mass relative to Sun. vx,vy,vz in AU/day." << endl;
  cout << "       Outputs total kinetic,potential energy, angular momentum," << endl;
  cout << "       x,y,z,vx,vy,vz,ep,ek,am, for each planet. vx,vy,vz in AU/yr." << endl;
  cout << "  -E   Euler's forward method with Earth in circular orbit aroudn the Sun" << endl;
  cout << "       Outputs x,y,z,vx,vy,vz, for each planet. vx,vy,vz in AU/yr." << endl;
  cout << "  -V   Velocity Verlet method with Earth in circular orbit aroudn the Sun" << endl;
  cout << "       Outputs x,y,z,vx,vy,vz, for each planet. vx,vy,vz in AU/yr." << endl;
  cout << "  -T   Time Euler's and Vel. Verlet with preset n=10^1 to 10^9. Outputs" << endl;
  cout << "       result averaged over 10 runs in csv format: n,euler,verlet" << endl;
  cout << "  -C   Creator mode. Experiment with Sun-Earth-Jupiter system. Select Earth's" << endl;
  cout << "       velocity [AU/yr.] with -v, include and specify Jupiter mass with -j," << endl;
  cout << "       play with the 1/r^b relationship between radius and force" << endl;
  cout << "       with -b." << endl;
  cout << "       Outputs total kinetic,potential energy, angular momentum," << endl;
  cout << "       x,y,z,vx,vy,vz,ep,ek,am, for each planet. vx,vy,vz in AU/yr." << endl;
  cout << "  -P   Calculates the perihelion precession of Mercury. Include gen. rel. with -r." << endl;
  cout << "  -n   Number of integration points per year." << endl;
  cout << "  -y   Years to simulate" << endl;
  cout << "  -b   1/r^b relationship between radius and gravitational force (only with -C)" << endl;
  cout << "  -v   Earth's velocity vy (only with -C)" << endl;
  cout << "  -j   Enables Jupiter under -C and specifies its mass multiplier" << endl;
  cout << "  -r   Includes general relativistic correction to Newton's grav. law. (only with -P)" << endl;
  cout << "" << endl;
  cout << "Example:" << endl;
  cout << "$ ./project3 -C -n 100 -y 2 -v 5" << endl;
  cout << "ek,ep,am,0x,0y,0z,0vx,0vy,0vz,0ep,0ek,0am,1x,1y,1z,1vx,1vy,1vz,1ep,1ek,1am," << endl;
  cout << "3.7500000000000003e-05,-1.1843525281307230e-04,1.5000000000000000e-05,0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,-1.1843525281307230e-04,0.0000000000000000e+00,0.0000000000000000e+00,1.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,5.0000000000000000e+00,0.0000000000000000e+00,-1.1843525281307230e-04,3.7500000000000003e-05,1.5000000000000000e-05" << endl;
  cout << "..." << endl;
}

/* From utils.cpp */


int main(int argc, char **argv) {
  int opt;
  int mode = -1;
  long long int n = 1;
  double years = 1.;
  double dt = 0.;
  double beta = 2.0;
  double earth_vy0 = 2*M_PI;
  double jupiter_mass = 954.79194E-6;  // relative to the sun
  bool include_jupiter = false;
  bool general_relativistic = false;  // enables gen. rel. grav. force in PERIHELION mode
  string fname = "";
  vec pos(3), vel(3);
  mat result, system;
  Solver solver;

  // Parsing args
  if (argc < 2) {
    cerr << "error: missing arguments" << endl;
    print_usage();
    exit(1);
   }

  opterr = 0;
  while ((opt = getopt(argc, argv, "hEVTCPn:y:f:b:v:j:r")) != -1) {
    switch (opt) {
    case 'h':
      print_usage();
      exit(0);
    case 'E':
      mode = EULER;
      continue;
    case 'V':
      mode = VERLET;
      continue;
    case 'T':
      mode = TIME;
      continue;
    case 'C':
      mode = CREATOR;
      continue;
    case 'P':
      mode = PERIHELION;
      continue;
    case 'n':
      n = (long long int) atoi(optarg);
      dt = 1./n;
      continue;
    case 'y':
      years = strtod(optarg, NULL);
      continue;
    case 'f':
      fname = optarg;
      mode = SYSTEM;
      continue;
    case 'b':
      beta = strtod(optarg, NULL);
      continue;
    case 'r':
      general_relativistic = true;
      continue;
    case 'v':
      earth_vy0 = strtod(optarg, NULL);
      continue;
    case 'j':
      jupiter_mass *= strtod(optarg, NULL);
      include_jupiter = true;
      continue;
    default:
      cerr << "error: unknown option: " << (char) optopt << endl;
      print_usage();
      exit(1);
    }
  }

  // Set number of steps to years times number of steps per year
  n = (long long int) n*years;

  // Run simulation
  switch (mode) {
    // Non-object oriented Euler
  case EULER:
    result = earth_circular_fwd_euler(dt, n);
    cout << "x,y,vx,vy" << endl;
    result.save(cout, csv_ascii);
    break;

    // Non-object oriented Verlet    
  case VERLET:
    result = earth_circular_verlet(dt, n);
    cout << "x,y,vx,vy" << endl;
    result.save(cout, csv_ascii);
    break;

    // Times Euler and Verlet with pre-defined time steps
  case TIME:
    result = time_algorithms();
    cout << "n,euler,verlet" << endl;
    result.save(cout, csv_ascii);
    break;

    // sun-earth system with fixed sun and variable grav. force, earth v0, size of Jupiter
  case CREATOR:
    {
      double ek0, ek, ep0, ep, am0, am;

      pos = {0.,0.,0.};
      vel = {0.,0.,0.};
      Planet sun = Planet(1., pos, vel);
      sun.fixed = true;
  
      pos = {1.,0.,0.};
      vel = {0., earth_vy0, 0.};
      Planet earth = Planet(3.0E-6, pos, vel);

      pos = {-5.2,0.,0.};
      vel = {0.,-2*M_PI/sqrt(5.2),0.};  // circular orbit
      Planet jupiter = Planet(jupiter_mass, pos, vel);

      solver.add(&sun);
      solver.add(&earth);
      if (include_jupiter)
	solver.add(&jupiter);

      am0 = earth.angular_momentum();
      ek0 = earth.kinetic_energy();
      ep0 = solver.potential_energy(&earth);

      solver.beta = beta;
      solver.solve(n, dt);

      am = earth.angular_momentum();
      ek = earth.kinetic_energy();
      ep = solver.potential_energy(&earth);

      cout << solver.csv_header() << endl;
      solver.flight_log.save(cout, csv_ascii);
    }
    break;

    // Calculate perihelion of planet Mercury
  case PERIHELION:
    cout << scientific;
    cout.precision(4);
    cout << perihelion_of_mercury(dt, n, general_relativistic) << endl;
    break;

    // Loads system from file and solves it
  case SYSTEM:
    system.load(fname, csv_ascii);
    solver.build(system);
    solver.solve(n, dt);
    cout << solver.csv_header() << endl;
    solver.flight_log.save(cout, csv_ascii);
    solver.nuke();
    break;

  default:
    break;
  }
  
  return 0;
}
