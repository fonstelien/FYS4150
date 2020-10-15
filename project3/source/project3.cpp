#include <fstream>
#include <unistd.h>
#include <time.h>
#include "utils.hpp"

/* Defaults */

/* Modes */

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
  vec euler_pos(3), euler_vel(3);
  vec verlet_pos(3), verlet_vel(3);
  int n = 0;
  double h = 0.;


  // Parsing args
  if (argc < 2) {
    cerr << "error: missing arguments" << endl;
    print_usage();
    exit(1);
   }

  opterr = 0;
  while ((opt = getopt(argc, argv, "hn:")) != -1) {
    switch (opt) {
    case 'h':
      print_usage();
      exit(0);
    case 'n':
      n = atoi(optarg);
      h = 1./n;
      continue;
    default:
      cerr << "error: unknown option: " << (char) optopt << endl;
      print_usage();
      exit(1);
    }
  }

  euler_pos[0] = 1;
  euler_pos[1] = 0;
  euler_pos[2] = 0;
  euler_vel[0] = 0;
  euler_vel[1] = 2*M_PI;
  euler_vel[2] = 0;

  verlet_pos[0] = 1;
  verlet_pos[1] = 0;
  verlet_pos[2] = 0;
  verlet_vel[0] = 0;
  verlet_vel[1] = 2*M_PI;
  verlet_vel[2] = 0;

  cout << "euler_pos0,euler_pos1,euler_pos2,verlet_pos0,verlet_pos1,verlet_pos2" << endl;
  for (int i = 0; i < n+1; i++) {
    cout << euler_pos[0] << "," << euler_pos[1] << "," << euler_pos[2] << ",";
    cout << verlet_pos[0] << "," << verlet_pos[1] << "," << verlet_pos[2] << endl;
    euler_fwd_iteration(h, euler_pos, euler_vel);
    velocity_verlet_iteration(h, verlet_pos, verlet_vel);
  }

  
  return 0;
}
