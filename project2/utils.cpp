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
