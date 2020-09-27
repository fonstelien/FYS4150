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

/* Helper function for generating general tridiagonal symmetric matrix of size n*n */
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

/* Finds the max fabs value along the off-diagonals in a symmetric matrix */
/* A and stores the position in *k, *l. Searches the lower diagonals. */
/* Returns fabs of max value. */
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

/* Finds and returns mat A's eigenvalues. The rotations are stored in mat S. */
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


class Polynomial {
protected:
  double h2;
  double d;
  double e;

public:
  virtual void init() {;}
  virtual double operator()(int k, double x) = 0;
};


class P: public Polynomial {
public:
  void init(double diag_element, double off_element) {
    d = diag_element;
    e = off_element;
  }

  double operator()(int k, double x) {
    int i;
    double pk_2, pk_1, pk;  // k minus 2, k minus 1, k
    
    pk_2 = d - x;
    if (k == 1)
      return pk_2;
    pk_1 = (d - x)*pk_2 - e*e;
    if (k == 2)
      return pk_1;

    pk = 0.;
    i = 2;
    while (i++ < k) {
      pk = (d - x)*pk_1 - e*e*pk_2;
      pk_2 = pk_1;
      pk_1 = pk;
    }
    
    return pk;
  }
};

class Q: public Polynomial {
public:
  void init(double diag_element, double off_element, double h) {
    h2 = h*h;
    d = diag_element;
    e = off_element;
  }

  double operator()(int k, double x) {
    int i;
    double qk_1, qk;  // k minus 1, k
    
    qk_1 = (d + h2) - x;
    if (k == 1)
      return qk_1;

    qk = 0.;
    i = 1;
    while (i++ < k) {
      qk = ((d + i*i*h2) - x) - e*e/qk_1;
      qk_1 = qk;
    }
    
    return qk;
  }
};


double bisection_root_finder(Polynomial &p, int k, double x_min, double x_max) {
  int i = 0;
  double x_left = x_min;
  double x_right = x_max;
  double x_mid;

  while ((x_right-x_left) > BISECTION_ROOT_FINDER_EPS &&	\
	 i++ < BISECTION_ROOT_FINDER_MAXITER) {

    x_mid = (x_right + x_left) / 2;
    if (p(k, x_left)*p(k, x_mid) < 0)
      x_right = x_mid;
    else
      x_left = x_mid;
  }

  return x_right;
}


vec quick_solver_tridiag_sym_toeplitz(int n) {
  double x_min, x_max, x_tmp;
  vec eigenvals(n);
  double h = 1./n;
  double d = 2.;
  double e = -1.;

  P p;
  p.init(d, e);

  x_max = d;
  x_min = 0.;
  eigenvals[0] = bisection_root_finder(p, 1, x_min, x_max);
  for (int k = 2; k < n+1; k++) {
    x_min = 0.;
    int i = 0;
    while (i < k-1) {
      x_max = eigenvals[i];
      eigenvals[i] = bisection_root_finder(p, k, x_min, x_max);
      x_min = x_max;
      i++;
    }
    x_min = x_max;
    x_max = d + fabs(2*e);
    eigenvals[i] = bisection_root_finder(p, k, x_min, x_max);
  }
  
  return eigenvals/(h*h);
}


vec quick_solver_tridiag_sym_general(int n, double rho_max) {
  double x_min, x_max, x_tmp;
  vec eigenvals(n);
  double h = rho_max/n;
  double h2 = h*h;
  double d = 2./h2;
  double e = -1./h2;
  Q q;

  q.init(d, e, h);
  
  x_max = d + h2;
  x_min = 0.;
  eigenvals[0] = bisection_root_finder(q, 1, x_min, x_max);
  for (int k = 2; k < n+1; k++) {
    x_min = 0.;
    int i = 0;
    while (i < k-1) {
      x_max = eigenvals[i];
      eigenvals[i] = bisection_root_finder(q, k, x_min, x_max);
      x_min = x_max;
      i++;
    }
    x_min = x_max;
    x_max = (d + k*k*h2) + fabs(2*e);
    eigenvals[i] = bisection_root_finder(q, k, x_min, x_max);
  }
  
  return eigenvals;
}
