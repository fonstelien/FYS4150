#include <iostream>
#include <armadillo>
#include <cmath>

#define DEBUG(msg) cout << "DEBUG " << msg << endl;

#define BISECTION_ROOT_FINDER_MAXITER 1E3
#define BISECTION_ROOT_FINDER_EPS 1.E-6


using namespace arma;
using namespace std;

/* From utils.cpp */
mat tridiag_sym_toeplitz(int n);
mat tridiag_sym_general(int n, double rho_max);
double max_off_diag_value(mat A, int *k, int *l);
vec jacobi_solver(mat A, mat &S, double tolerance);
void sort_eigen(vec &eigenvals, mat &eigenvecs);
vec quick_solver(int n, double rho_max);


