#include <iostream>
#include <armadillo>
#include <cmath>
#include <unistd.h>
#include <time.h>

#define DEBUG(msg) cout << "DEBUG " << msg << endl;

using namespace arma;
using namespace std;

/* From utils.cpp */
mat tridiag_sym_toeplitz(int n, double center_element, double off_element);
double max_off_diag_value(mat A, int *k, int *l);
vec jacobi_solver(mat A, mat &S, double tolerance);
void sort_eigen(vec &eigenvals, mat &eigenvecs);


