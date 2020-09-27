#include <iostream>
#include <armadillo>
#include <cmath>

#define DEBUG(msg) cout << "DEBUG " << msg << endl;

#define BISECTION_ROOT_FINDER_MAXITER 50
#define BISECTION_ROOT_FINDER_EPS 1.E-3

using namespace arma;
using namespace std;

/* From utils.cpp */
mat make_tridiag_sym_toeplitz(int n);
mat make_tridiag_sym_general(int n, double rho_max);
double max_off_diag_value(mat A, int *k, int *l);
vec jacobi_solver(mat A, mat &S, double tolerance);
void sort_eigen_pairs(vec &eigenvals, mat &eigenvecs);
vec quick_solver_tridiag_sym_toeplitz(int n);
vec quick_solver_tridiag_sym_general(int n, double rho_max);
