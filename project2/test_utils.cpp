#include "test_utils.hpp"

using namespace arma;
using namespace std;

vec tridiag_sym_toeplitz_exact_eigenvals(mat A) {
  int n = A.n_rows;
  vec eigenvals = vec(n);
  double d = A(0,0);
  double a = A(0,1);

  for (int i = 0; i < n; i++)
    eigenvals[i] = d + 2*a*cos((i+1)*M_PI/(n+1));

  return eigenvals;
}

mat tridiag_sym_toeplitz_exact_eigenvecs(mat A) {
  int n = A.n_rows;
  mat eigenvecs = mat(arma::size(A));

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      eigenvecs(i,j) = sin((i+1)*(j+1)*M_PI/(n+1));

  return eigenvecs;
}
