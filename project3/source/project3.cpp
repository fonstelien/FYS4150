#include <fstream>
#include <unistd.h>
#include <time.h>
#include "utils.hpp"

/* Defaults */

/* Modes */
#define EULER 1
#define VERLET 2
#define SYSTEM 3

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
  string fname;
  double dt = 0.;
  vec pos(3), vel(3);
  mat flight_log, system;
  Solver solver;

  // Parsing args
  if (argc < 2) {
    cerr << "error: missing arguments" << endl;
    print_usage();
    exit(1);
   }

  opterr = 0;
  while ((opt = getopt(argc, argv, "hevn:y:f:")) != -1) {
    switch (opt) {
    case 'h':
      print_usage();
      exit(0);
    case 'e':
      mode = EULER;
      continue;
    case 'v':
      mode = VERLET;
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
    flight_log = mat(n, 6);
    earth_circular_fwd_euler(dt, flight_log);
    cout << "x,y,z,vx,vy,vz" << endl;
    flight_log.save(cout, csv_ascii);
    break;
    
  case VERLET:
    flight_log = mat(n, 6);    
    earth_circular_verlet(dt, flight_log);
    cout << "x,y,z,vx,vy,vz" << endl;
    flight_log.save(cout, csv_ascii);
    break;

  case SYSTEM:
    system.load(fname, csv_ascii);
    solver.build(system);
    solver.solve(n, dt);
    cout << solver.csv_header() << endl;
    solver.flight_log.save(cout, csv_ascii);
    break;

  default:
    break;
  }
  
  return 0;
}
