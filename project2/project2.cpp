#include "utils.hpp"

#define TOL 1.E-9

using namespace arma;
using namespace std;

/* From utils.cpp */
mat tridiag_sym_toeplitz(int n, double center_element, double off_element);
double max_off_diag_value(mat A, int *k, int *l);
vec jacobi_solver(mat A, mat &S, double tolerance);
void sort_eigen(vec &eigenvals, mat &eigenvecs);


int main() {
  int n = 10;
  mat A = tridiag_sym_toeplitz(n, 2, -1);
  mat S(n, n, fill::eye);
  vec eigenvals;
  
  eigenvals = jacobi_solver(A, S, TOL);

  cout << eigenvals << endl;
  cout << S << endl;

  sort_eigen(eigenvals, S);

  cout << eigenvals << endl;
  cout << S << endl;
  
  return 0;
}
