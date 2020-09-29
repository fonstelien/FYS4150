#include "utils.hpp"

using namespace arma;
using namespace std;


/* Helper function for generating tridiagonal symmetric Toeplitz matrix of size n*n */
mat make_tridiag_sym_toeplitz(int n) {
  double h = 1./n;
  double center_element = 2/(h*h);
  double off_element = -1/(h*h);
  mat A = zeros<mat>(n,n);
  A.diag(0) = zeros<rowvec>(n) + center_element;
  A.diag(1) = zeros<rowvec>(n-1) + off_element;
  A.diag(-1) = zeros<rowvec>(n-1) + off_element;

  return A;
}


/* Helper function for generating general tridiagonal symmetric matrix of size n*n, */
/* With constant values along the off-diagonals and row-dependent value along the central. */
mat make_tridiag_sym_general(int n, double rho_max) {
  double h = rho_max/n;
  double off_element = -1/(h*h);
  mat A = zeros<mat>(n,n);

  for (int i = 0; i < n; i++)
    A(i,i) = 2/(h*h) + (i*h)*(i*h);

  A.diag(1) = zeros<rowvec>(n-1) + off_element;
  A.diag(-1) = zeros<rowvec>(n-1) + off_element;

  return A;
}


/* Finds the max fabs() value along the off-diagonals in a symmetric matrix */
/* A and stores the position in *k, *l. Searches only the lower diagonals. */
/* Returns fabs() of max value. */
double max_off_diag_value(mat A, int *k, int *l) {
  int n = A.n_rows;
  double max_val, new_max;

  *k = 1;
  *l = 0;
  max_val = fabs(A(*k,*l));

  for (int i = 2; i < n; i++)
    for (int j = 0; j < i; j++)
      if ((new_max = fabs(A(i,j))) > max_val) {
	max_val = new_max;
	*k = i;
	*l = j;
      }

  return max_val;
}


/* Finds and returns mat A's eigenvalues by the Jacobi method. */
/* The rotations are stored in mat S. */
vec jacobi_solver(mat A, mat &S, double tolerance) {
  double a_ik, a_il, a_kk, a_kl, a_ll;
  double s_ki, s_li;
  double tau, t, s, c;  // t=tan(theta), s=sin(theta), c=cos(theta)
  int k, l;
  int n = A.n_rows;
  vec eigenvals(n);

  /* Rotating */
  while (max_off_diag_value(A, &k, &l) > tolerance) {  
    /* Finding the right rotation */
    a_kk = A(k,k);
    a_kl = A(k,l);
    a_ll = A(l,l);
  
    tau = (a_ll - a_kk) / (2*a_kl);
    if (tau > 0)
      t = 1/(tau + sqrt(tau*tau + 1));
    else
      t = 1/(tau - sqrt(tau*tau + 1));

    c = 1/sqrt(1 + t*t);
    s = c*t;

    /* Updating A's values*/
    A(k,k) = a_kk*c*c - 2*a_kl*s*c + a_ll*s*s;
    A(k,l) = A(l,k) = 0.;
    A(l,l) = a_ll*c*c + 2*a_kl*s*c + a_kk*s*s;

    for (int i = 0; i < n; i++) {
      if (i == k || i == l)
	continue;
    
      a_ik = A(i,k);
      a_il = A(i,l);
      A(i,k) = A(k,i) = a_ik*c - a_il*s;
      A(i,l) = A(l,i) = a_il*c + a_ik*s;
    }

    /* Updating S's values */
    for (int i = 0; i < n; i++) {
      s_ki = S(k,i);
      s_li = S(l,i);
      S(k,i) = s_ki*c - s_li*s;
      S(l,i) = s_ki*s + s_li*c;
    }
  }

  /* Copy eigenvalues from A into eigenvals */
  for (int i = 0; i < n; i++)
    eigenvals[i] = A(i,i);

  return eigenvals;
}
