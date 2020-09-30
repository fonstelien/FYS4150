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
#define PLOT_VECTORS 3
#define TIME 4

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
  cout << "usage: project2 [-h | -C | -I | -V | -T] [-j | -p] [-t | -g] [-re] n" << endl;
  cout << "  -h   print this text" << endl;
  cout << "  -C   Compare exact and numerical results. Prints to stdout." << endl;
  cout << "  -I   number of Iterations before convergence (Jacobi only)" << endl;
  cout << "  -V   prints eigenVectors to stdout (Jacobi only)" << endl;
  cout << "  -T   prints CPU Time for algorithm" << endl;
  cout << "  -j   run Jacobi algorithm" << endl;
  cout << "  -p   run Polynomial expansion algorithm" << endl;
  cout << "  -t   Toeplitz tridiagonal matrix" << endl;
  cout << "  -g   'general' tridiagonal matrix" << endl;
  cout << "  -r   rho max in the 'general' matrix" << endl;
  cout << "  -e   convergence tolerance (eps)" << endl;
  cout << "  -n   " << endl;
  cout << "Where suitable, the results are printed in csv format.";
  cout << "example:";
  cout << "$ project2 -Cjg -e 1.E-2 -r 10.0 -n 100" << endl;
  cout << "exact,numeric" << endl;
  cout << "3.000000000000e+00,2.779213340881e+00" << endl;
  cout << "7.000000000000e+00,6.654696896637e+00" << endl;
  cout << "1.100000000000e+01,1.054835747359e+01" << endl;
  cout << "..." << endl;
}

/* From utils.cpp */
void sort_eigen_pairs(vec &eigenvals, mat &eigenvecs);
void sort_eigenvals(vec &eigenvals);
vec tridiag_sym_toeplitz_exact_eigenvals(int n);
vec tridiag_sym_general_exact_eigenvals(int n);
mat tridiag_sym_toeplitz_exact_eigenvecs(int n);

/* From jacobi.cpp */
mat make_tridiag_sym_toeplitz(int n);
mat make_tridiag_sym_general(int n, double rho_max);
double max_off_diag_value(mat A, int *k, int *l);
vec jacobi_solver(mat A, mat &S, int *rotations, double tolerance);
  
/* From polynomial_expansion.cpp */
vec poly_exp_solver_tridiag_sym_toeplitz(int n, double eps);
vec poly_exp_solver_tridiag_sym_general(int n, double rho_max, double eps);


int main(int argc, char **argv) {
  vec eigenvals_numeric, eigenvals_exact;
  mat A, S, eigenvecs_exact;
  int n = DIMENSION;
  float h;
  float multiplier = 0;
  ofstream out_file;
  clock_t start_t = 0;
  clock_t end_t = 0;
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
  while ((opt = getopt(argc, argv, "hCIVTjptgr:e:n:")) != -1) {
    switch (opt) {
    case 'h':
      print_usage();
      exit(0);
    case 'n':
      n = atoi(optarg);
      continue;
    case 'C':
      mode = COMPARE;
      continue;
    case 'I':
      mode = ITERATIONS;
      continue;
    case 'V':
      mode = PLOT_VECTORS;
      algorithm = JACOBI;
      matrix = TOEPLITZ;
      continue;
    case 'T':
      mode = TIME;
      continue;
    case 'j':
      algorithm = JACOBI;
      continue;
    case 'p':
      algorithm = POLY_EXP;
      continue;
    case 't':
      matrix = TOEPLITZ;
      continue;
    case 'g':
      matrix = GENERAL;
      continue;
    case 'r':
      rho_max = strtod(optarg, NULL);
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
    start_t = clock();
    eigenvals_numeric = jacobi_solver(A, S, &iterations, eps);
    end_t = clock();
    sort_eigen_pairs(eigenvals_numeric, S);
    break;

  case POLY_EXP:
    start_t = clock();
    if (matrix == TOEPLITZ)
      eigenvals_numeric = poly_exp_solver_tridiag_sym_toeplitz(n, eps);
    else if (matrix == GENERAL)
      eigenvals_numeric = poly_exp_solver_tridiag_sym_general(n, rho_max, eps);
    end_t = clock();
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
  else if (mode == PLOT_VECTORS) {
    eigenvals_exact = tridiag_sym_toeplitz_exact_eigenvals(n);
    eigenvecs_exact = tridiag_sym_toeplitz_exact_eigenvecs(n);
    sort_eigen_pairs(eigenvals_exact, eigenvecs_exact);
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
    cout << n << "," << iterations << endl;
    break;
  case PLOT_VECTORS:
    {
      vec u = eigenvecs_exact.col(0);
      mat V = S * eigenvecs_exact;
      vec v = V.col(0);
      cout.precision(3);
      cout << scientific;
      cout << "exact,rotated";
      for (int i = 0; i < n; i++)
	cout << ",s" << i;
      cout << endl;
      for (int i = 0; i < n; i++) {
      	cout << u[i] << "," << v[i];
	for (int j = 0; j < n; j++)
	  cout << "," << S(i,j);
	cout << endl;
      }
    }
    break;
  case TIME:
    lapsed_time = (double) (end_t - start_t) / CLOCKS_PER_SEC;
    cout << lapsed_time << endl;
  default:
    break;
  }

  
  return 0;
}
