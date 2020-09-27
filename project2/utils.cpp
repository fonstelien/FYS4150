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


/* Abstract class */
class Polynomial {
protected:
  double h2;
  double d;
  double e;

public:
  virtual void init() {;}
  virtual double operator()(int k, double x) = 0;
};


/* Class for solving the eigenvalue problem by expanding the polynomial */
/* P2(x) = (d-x)*P1(x) - e^2*P0(x) */
class P: public Polynomial {
public:
  void init(double d0, double e0) {
    d = d0;
    e = e0;
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


/* Class for solving the eigenvalue problem by expanding the polynomial */
/* Q2(x) = (d-x) - e^2/Q1(x) */
class Q: public Polynomial {
  double di(int i) {
    return d + i*i*h2;
  }
  
public:
  void init(double d0, double e0, double h) {
    h2 = h*h;
    d = d0;
    e = e0;
  }

  double operator()(int k, double x) {
    int i;
    double qk_1, qk;  // k minus 1, k

    i = 1;
    qk = di(i) - x;
    if (k == 1)
      return qk;

    qk_1 = qk;
    while (i++ < k) {
      qk = (di(i) - x) - e*e/qk_1;
      qk_1 = qk;
    }
    
    return qk;
  }
};


/* Finds the root of Pk(x) in the range [x_min, x_max] by the bisection method. */
double bisection_root_finder(Polynomial &p, int k, double x_min, double x_max) {
  int i = 0;
  double x_left = x_min;
  double x_right = x_max;
  double x_mid;
  double p_mid;
  
  while ((x_right-x_left)/(x_max-x_min) > BISECTION_ROOT_FINDER_EPS &&	\
	 i++ < BISECTION_ROOT_FINDER_MAXITER) {
    x_mid = (x_right + x_left) / 2;
    p_mid = p(k, x_mid);
    
    if (p_mid == 0)
      return x_mid;

    if (p(k, x_left)*p_mid < 0)
      x_right = x_mid;
    else
      x_left = x_mid;
  }

  return x_right;
}


/* Quick solver for tridiagonal symmetrical Toeplitz matrix. Finds the roots by polynomial */
/* expansion and bisection root search. */
vec quick_solver_tridiag_sym_toeplitz(int n) {
  double x_min, x_max;
  vec eigenvals(n);
  double h = 1./n;
  double d = 2.;
  double e = -1.;

  /* P2(x) = (d-x)*P1(x) - e^2*P0(x) */  
  P p;
  p.init(d, e);

  /* Finding roots by polynomial expansion */  
  eigenvals[0] = d - sqrt(-e);  // roots for k=2
  eigenvals[1] = d + sqrt(-e);
  for (int k = 3; k < n+1; k++) {
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


/* Quick solver for tridiagonal symmetrical Toeplitz matrix. Finds the roots by polynomial */
/* expansion and bisection root search. */
vec quick_solver_tridiag_sym_general(int n, double rho_max) {
  double x_min, x_max;
  vec eigenvals(n);
  double h = rho_max/n;
  double h2 = h*h;
  double d = 2./h2;
  double e = -1./h2;

  /* Q2(x) = (d-x) - e^2/Q1(x) */
  Q q;
  q.init(d, e, h);

  /* Finding roots by polynomial expansion */  
  eigenvals[0] = d + h2;  // root for k=1
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
