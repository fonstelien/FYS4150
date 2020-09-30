#include "catch.hpp"
#include "utils.hpp"

using namespace arma;
using namespace std;

#define TOL 1.E-6

/* From utils.cpp */
void sort_eigen_pairs(vec &eigenvals, mat &eigenvecs);
void sort_eigenvals(vec &eigenvals);
vec tridiag_sym_toeplitz_exact_eigenvals(int n);
vec tridiag_sym_general_exact_eigenvals(int n);
mat tridiag_sym_toeplitz_exact_eigenvecs(int n);

/* From jacobi.cpp */
mat make_tridiag_sym_toeplitz(int n);
mat make_tridiag_sym_general(int n, double rho_max);
double max_off_diag_value(mat A, int *k, int *l);
vec jacobi_solver(mat A, mat &S, int *rotations, double tolerance);
  
/* From polynomial_expansion.cpp */
vec poly_exp_solver_tridiag_sym_toeplitz(int n, double eps);
vec poly_exp_solver_tridiag_sym_general(int n, double rho_max, double eps);


TEST_CASE("Max value on off-diagonal", "[max-off-diag]") {
  double set_val, max_val;
  int i,j,k,l;
  int n = 10;
  mat A = zeros<mat>(n,n);

  i = 4;
  j = 3;
  set_val = 5.;
  A(i,j) = set_val;
  max_val = max_off_diag_value(A, &k, &l);
  REQUIRE(i == k);
  REQUIRE(j == l);
  REQUIRE(max_val == fabs(set_val));

  i = 7;
  j = 2;
  set_val = -6.;
  A(i,j) = set_val;
  max_val = max_off_diag_value(A, &k, &l);
  REQUIRE(i == k);
  REQUIRE(j == l);
  REQUIRE(max_val == fabs(set_val));
}


TEST_CASE("Toeplitz tridiag eigen pairs", "[toeplitz-eigen-pairs]") {
  vec eigenvals_numerical, eigenvals_exact;
  mat eigenvecs_numerical, eigenvecs_exact;
  int n, rotations;
  mat A, S, U, V;

  n = 10;
  
  // eig_sym(eigenvals_numerical, eigenvecs_numerical, A);
  eigenvals_exact = tridiag_sym_toeplitz_exact_eigenvals(n);
  eigenvecs_exact = tridiag_sym_toeplitz_exact_eigenvecs(n);
  sort_eigen_pairs(eigenvals_exact, eigenvecs_exact);
  V = eigenvecs_exact;

  A = make_tridiag_sym_toeplitz(n);
  S.eye(n,n);  
  eigenvals_numerical = jacobi_solver(A, S, &rotations, 1.E-9);
  sort_eigen_pairs(eigenvals_numerical, S);
  U = S * V;
    
  for (int i = 0; i < n; i++)
    REQUIRE(eigenvals_numerical[i] == Approx(eigenvals_exact[i]).epsilon(TOL));

  for (int i = 0; i < n; i++) {
    vec u_i = U.col(i);
    vec v_i = V.col(i);
    for (int j = 0; j < n; j++) {
      vec u_j = U.col(j);
      vec v_j = V.col(j);
      REQUIRE(fabs(dot(u_i, u_j)) == Approx(fabs(dot(v_i, v_j))).margin(TOL));
    }
  }
}

TEST_CASE("Quick Toeplitz", "[quick-toeplitz]") {
  vec eigenvals_numerical, eigenvals_exact;
  int n;
  mat A;

  n = 100;
  eigenvals_exact = tridiag_sym_toeplitz_exact_eigenvals(n);
  sort_eigenvals(eigenvals_exact);

  eigenvals_numerical = poly_exp_solver_tridiag_sym_toeplitz(n, 1.E-6);
  
  for (int i = 0; i < n; i++)
    REQUIRE(eigenvals_numerical[i] == Approx(eigenvals_exact[i]).epsilon(1.E-6));
}

TEST_CASE("Quick general", "[quick-general]") {
  vec eigenvals_numerical, eigenvals_exact;
  int N, n;
  mat A;
  double rho_max;
  
  N = 3;
  eigenvals_exact = tridiag_sym_general_exact_eigenvals(N);

  n = 500;    
  rho_max = 4.;
  eigenvals_numerical = poly_exp_solver_tridiag_sym_general(n, rho_max, 1.E-6);
  
  for (int i = 0; i < N; i++)
    REQUIRE(eigenvals_numerical[i] == Approx(eigenvals_exact[i]).epsilon(1.E-2));
}
