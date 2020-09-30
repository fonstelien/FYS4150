#include "utils.hpp"

using namespace arma;
using namespace std;

#define BISECTION_ROOT_FINDER_MAXITER 50

/* Global variable! */
double bisection_root_finder_eps;


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
  
  while ((x_right-x_left)/(x_max-x_min) > bisection_root_finder_eps &&	\
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
vec poly_exp_solver_tridiag_sym_toeplitz(int n, double eps) {
  double x_min, x_max;
  vec eigenvals(n);
  double h = 1./n;
  double d = 2.;
  double e = -1.;

  /* Global variable! */
  bisection_root_finder_eps = eps;

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
vec poly_exp_solver_tridiag_sym_general(int n, double rho_max, double eps) {
  double x_min, x_max;
  vec eigenvals(n);
  double h = rho_max/n;
  double h2 = h*h;
  double d = 2./h2;
  double e = -1./h2;

  /* Global variable! */
  bisection_root_finder_eps = eps;

  
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
