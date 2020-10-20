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
  cout << "usage: project3" << endl;
  cout << "... pending ..." << endl;
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
  bool stability = false;  // preservation of Ep+Ek, angular momentum
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
  while ((opt = getopt(argc, argv, "hEVTCPsrn:y:f:b:v:j:")) != -1) {
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
    case 's':
      stability = true;
      continue;
    case 'n':
      n = atoi(optarg);
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

  // Set number of steps to years times number of steps per year, +1 for the initial state
  n = (int) n*years + 1;

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
      ep0 = solver.potential_energy(&earth, &sun);
      if (include_jupiter)
	ep0 += solver.potential_energy(&earth, &jupiter);

      solver.beta = beta;
      solver.solve(n, dt);

      am = earth.angular_momentum();
      ek = earth.kinetic_energy();
      ep = solver.potential_energy(&earth, &sun);
      if (include_jupiter)
	ep += solver.potential_energy(&earth, &jupiter);

      if (stability) {
	cout << "ek0+ep0,ek+ep,am0,am" << endl;
	cout << ek0+ep0 << "," << ek+ep << "," << am0 << "," << am << endl;
      }
      else {
	cout << solver.csv_header() << endl;
	solver.flight_log.save(cout, csv_ascii);
      }
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
