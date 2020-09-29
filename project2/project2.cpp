#include <fstream>
#include <unistd.h>
#include <time.h>
#include "utils.hpp"

/* Defaults */
#define EPS 1.E-6  // tolerance
#define DIMENSION 10  // matrix dimension

/* Modes */
#define COMPARE 1
#define ITERATIONS 2
#define TIME 3

/* Algorithms */
#define JACOBI 1
#define POLY_EXP 2

/* Toeplitz or general matrix */
#define TOEPLITZ 1
#define GENERAL 2


using namespace arma;
using namespace std;


/* Prints program usage */
void print_usage() {
  cout << "usage: project1 [-h | -c | -t] [-jpe] n" << endl;
  cout << "  -h   print this text" << endl;
}


/* From utils.cpp */


int main(int argc, char **argv) {
  vec eigenvals_numeric, eigenvals_exact;
  mat A, S;
  int n = DIMENSION;
  float h;
  float multiplier = 0;
  ofstream out_file;
  clock_t start_t, end_t;
  double lapsed_time;
  int opt;
  int mode = -1;
  int algorithm = -1;
  int matrix = -1;
  double rho_max = 1.;
  double eps = EPS;
  int iterations;
  
  // Parsing args
  if (argc < 2) {
    cerr << "error: missing arguments" << endl;
    print_usage();
    exit(1);
   }

  opterr = 0;
  while ((opt = getopt(argc, argv, "hn:cjtgie:")) != -1) {
    switch (opt) {
    case 'h':
      print_usage();
      exit(0);
    case 'n':
      n = atoi(optarg);
      continue;
    case 'c':
      mode = COMPARE;
      continue;
    case 'i':
      mode = ITERATIONS;
      continue;
    case 'j':
      algorithm = JACOBI;
      continue;
    case 't':
      matrix = TOEPLITZ;
      continue;
    case 'g':
      matrix = GENERAL;
      continue;
    case 'e':
      eps = strtod(optarg, NULL);
      continue;
    default:
      cerr << "error: unknown option: " << (char) optopt << endl;
      print_usage();
      exit(1);
    }
  }
  
  h = rho_max/n;

  /* Running algorithms */
  switch (algorithm) {
  case JACOBI:
    if (matrix == TOEPLITZ)
      A = make_tridiag_sym_toeplitz(n);
    else if (matrix == GENERAL)
      A = make_tridiag_sym_general(n, rho_max);
    S.eye(n,n);
    eigenvals_numeric = jacobi_solver(A, S, &iterations, eps);
    sort_eigen_pairs(eigenvals_numeric, S);
    break;

  case POLY_EXP:
    if (matrix == TOEPLITZ)
      ;
    else if (matrix == GENERAL)
      ;
    break;
    
  default:
    break;
  }
  
  if (mode == COMPARE) {
    if (matrix == TOEPLITZ)
      eigenvals_exact = tridiag_sym_toeplitz_exact_eigenvals(n);
    else if (matrix == GENERAL)
      eigenvals_exact = tridiag_sym_general_exact_eigenvals(n);
    sort_eigenvals(eigenvals_exact);
  }

  /* Print results */
  switch (mode) {
  case COMPARE:
    cout.precision(12);
    cout << scientific;
    cout << "exact" << "," << "numeric" << endl;
    for (int i = 0; i < n; i++)
      cout << eigenvals_exact[i] << "," << eigenvals_numeric[i] << endl;
    break;

  case ITERATIONS:
    cout << iterations << endl;

  default:
    break;
  }

  
  return 0;
}
