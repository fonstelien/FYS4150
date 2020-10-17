#include <fstream>
#include <unistd.h>
#include <time.h>
#include "utils.hpp"

/* Defaults */

/* Modes */
#define EULER 1
#define VERLET 2
#define TIME 3
#define FORCE 4
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
  int n = 1;
  double years = 1.;
  double dt = 0.;
  double beta = 2.;
  string fname;
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
  while ((opt = getopt(argc, argv, "hEVTFn:y:f:b:")) != -1) {
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
    case 'F':
      mode = FORCE;
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
  case EULER:
    result = earth_circular_fwd_euler(dt, n);
    cout << "x,y,z,vx,vy,vz" << endl;
    result.save(cout, csv_ascii);
    break;
    
  case VERLET:
    result = earth_circular_verlet(dt, n);
    cout << "x,y,z,vx,vy,vz" << endl;
    result.save(cout, csv_ascii);
    break;

  case TIME:
    result = time_algorithms();
    cout << "n,euler,verlet" << endl;
    result.save(cout, csv_ascii);
    break;

  case FORCE:
    {
      pos = {0.,0.,0.};
      vel = {0.,0.,0.};
      Planet sun = Planet(1., pos, vel);
      sun.fixed = true;
  
      pos = {1.,0.,0.};
      vel = {0.,5.,0.};  // elliptic orbit
      Planet earth = Planet(3.0E-6, pos, vel);

      solver.beta = beta;
      solver.add(&sun);
      solver.add(&earth);
      solver.solve(n, dt);
      result = solver.flight_log.cols(6,11);
      cout << "x,y,z,vx,vy,vz" << endl;
      result.save(cout, csv_ascii);
    }
    break;

    
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
