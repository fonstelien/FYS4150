#include "utils.hpp"

using namespace arma;
using namespace std;


/* Sorts eigenpairs in ascending order with respect to eigenvalues */
void sort_eigen_pairs(vec &eigenvals, mat &eigenvecs) {
  vec tmp_vec;
  double tmp_val;
  int n = eigenvecs.n_cols;

  /* Sort ascending eigenvecs and eigenvals synchronously */
  int i, j;
  i = j = 1;
  while (i < n) {
    j = i;
    while (j > 0 && eigenvals[j] < eigenvals[j-1]) {
      // swap vals
      tmp_val = eigenvals[j-1];
      eigenvals[j-1] = eigenvals[j];
      eigenvals[j] = tmp_val;

      // swap vecs
      tmp_vec = eigenvecs.col(j-1);
      eigenvecs.col(j-1) = eigenvecs.col(j);
      eigenvecs.col(j) = tmp_vec;
      j--;
    }
    i++;
  }
}


/* Sorts eigenpairs in ascending order with respect to eigenvalues */
void sort_eigenvals(vec &eigenvals) {
  double tmp_val;
  int n = eigenvals.n_elem;

  /* Sort ascending eigenvecs and eigenvals synchronously */
  int i, j;
  i = j = 1;
  while (i < n) {
    j = i;
    while (j > 0 && eigenvals[j] < eigenvals[j-1]) {
      tmp_val = eigenvals[j-1];
      eigenvals[j-1] = eigenvals[j];
      eigenvals[j] = tmp_val;
      j--;
    }
    i++;
  }
}



vec tridiag_sym_toeplitz_exact_eigenvals(int n) {
  vec eigenvals = vec(n);
  double h = 1./n;
  double d = 2./(h*h);
  double a = -1./(h*h);

  for (int i = 0; i < n; i++)
    eigenvals[i] = d + 2*a*cos((i+1)*M_PI/(n+1));

  return eigenvals;
}

mat tridiag_sym_toeplitz_exact_eigenvecs(int n) {
  mat eigenvecs(n,n);

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      eigenvecs(i,j) = sin((i+1)*(j+1)*M_PI/(n+1));

  return eigenvecs;
}

vec tridiag_sym_general_exact_eigenvals(int n) {
  vec eigenvals = vec(n);

  eigenvals[0] = 3.;
  for (int i = 1; i < n; i++)
    eigenvals[i] = eigenvals[i-1] + 4.;

  return eigenvals; 
}
