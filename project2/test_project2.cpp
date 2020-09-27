#include "catch.hpp"
#include "test_utils.hpp"

using namespace arma;
using namespace std;

#define TOL 1.E-6

/* From test_utils.cpp */

/* From utils.cpp */



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
  int n;
  mat A, S, U, V;

  n = 10;
  A = make_tridiag_sym_toeplitz(n);
  S.eye(n,n);
  
  // eig_sym(eigenvals_numerical, eigenvecs_numerical, A);
  eigenvals_exact = tridiag_sym_toeplitz_exact_eigenvals(A);
  eigenvecs_exact = tridiag_sym_toeplitz_exact_eigenvecs(A);
  sort_eigen_pairs(eigenvals_exact, eigenvecs_exact);
  V = eigenvecs_exact;

  eigenvals_numerical = jacobi_solver(A, S, 1.E-9);
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
  A = make_tridiag_sym_toeplitz(n);
  eigenvals_exact = tridiag_sym_toeplitz_exact_eigenvals(A);
  sort_eigen_pairs(eigenvals_exact, A);

  eigenvals_numerical = quick_solver_tridiag_sym_toeplitz(n);
  
  for (int i = 0; i < n; i++)
    REQUIRE(eigenvals_numerical[i] == Approx(eigenvals_exact[i]).epsilon(1.E-6));
}

TEST_CASE("Quick general", "[quick-general]") {
  vec eigenvals_numerical, eigenvals_exact;
  int n;
  mat A;
  double rho_max;
  
  n = 500;
  eigenvals_exact = vec(3);  
  eigenvals_exact[0] = 3.;
  for (int i = 1; i < 3; i++)
    eigenvals_exact[i] = eigenvals_exact[i-1] + 4.;

  rho_max = 4.;
  eigenvals_numerical = quick_solver_tridiag_sym_general(n, rho_max);
  
  for (int i = 0; i < 3; i++)
    REQUIRE(eigenvals_numerical[i] == Approx(eigenvals_exact[i]).epsilon(1.E-2));
}
