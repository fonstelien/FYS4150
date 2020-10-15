#include <fstream>
#include <unistd.h>
#include <time.h>
#include "utils.hpp"

/* Defaults */

/* Modes */
#define EULER 1
#define VERLET 2

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
  double h = 0.;
  double years = 1.;
  vec pos(3), vel(3);
  mat results;

  // Parsing args
  if (argc < 2) {
    cerr << "error: missing arguments" << endl;
    print_usage();
    exit(1);
   }

  opterr = 0;
  while ((opt = getopt(argc, argv, "hevn:y:")) != -1) {
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
      h = 1./n;
      continue;
    case 'y':
      years = strtod(optarg, NULL);
      continue;
    default:
      cerr << "error: unknown option: " << (char) optopt << endl;
      print_usage();
      exit(1);
    }
  }

  /* Set number of steps to years times number of steps per year, +1 for the initial state */
  n = (int) n*years + 1;

  pos = {0.,0.,0.};
  vel = {0.,0.,0.};
  Planet sun = Planet(1, pos, vel);
  pos = {1.,0.,0.};
  vel = {0.,2*M_PI,0.};  
  Planet earth = Planet(6.0E24/2.0E30, pos, vel);

  pos = {2.550209154629485E+00,-4.432721232593654E+00,-3.866707508925721E-02};
  vel = {6.447098824156304E-03,4.121019457101368E-03,-1.613529989591600E-04};
  vel *= 365;
  Planet jupiter = Planet(1.9E27/2.0E30, pos, vel);

  
  Solver solver;
  solver.add(sun);
  solver.add(earth);
  solver.add(jupiter);
  solver.solve(n, h);
  solver.to_csv();
  
  if (mode == EULER || mode == VERLET)
    results = mat(n, 6);  // postitions and velocities for all n
 
  switch (mode) {
  case EULER:
    earth_circular_fwd_euler(h, results);
    cout << "x,y,z,vx,vy,vz" << endl;
    results.save(cout, csv_ascii);
    break;
  case VERLET:
    earth_circular_verlet(h, results);
    cout << "x,y,z,vx,vy,vz" << endl;
    results.save(cout, csv_ascii);
    break;
  default:
    break;
  }

  
  return 0;
}
