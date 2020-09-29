#include <iostream>
#include <armadillo>
#include <cmath>

#define DEBUG(msg) cout << "DEBUG " << msg << endl;

#define BISECTION_ROOT_FINDER_MAXITER 50
#define BISECTION_ROOT_FINDER_EPS 1.E-6

using namespace arma;
using namespace std;

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
vec quick_solver_tridiag_sym_toeplitz(int n);
vec quick_solver_tridiag_sym_general(int n, double rho_max);
