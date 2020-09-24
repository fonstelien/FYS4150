#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "test_utils.hpp"

using namespace arma;
using namespace std;

#define TOL 1.E-6

/* From test_utils.cpp */
vec tridiag_sym_toeplitz_exact_eigenvals(mat A);
mat tridiag_sym_toeplitz_exact_eigenvecs(mat A);

/* From utils.cpp */
mat tridiag_sym_toeplitz(int n, double center_element, double off_element);
double max_off_diag_value(mat A, int *k, int *l);
vec jacobi_solver(mat A, mat &S, double tolerance);
void sort_eigen(vec &eigenvals, mat &eigenvecs);


TEST_CASE("Sym. tridiag. max value on off-diagonal", "[sym-tri-max]") {
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




TEST_CASE("Sym. tridiag. Toeplitz eigen", "[sym-tri-toeplitz-eigen]") {
  vec eigenvals_numerical, eigenvals_exact;
  mat eigenvecs_numerical, eigenvecs_exact;
  int n = 5;
  mat A = tridiag_sym_toeplitz(n, 2., -1.);
  mat S(n, n, fill::eye);
  mat U1, U2;
  
  // eig_sym(eigenvals_numerical, eigenvecs_numerical, A);
  eigenvals_exact = tridiag_sym_toeplitz_exact_eigenvals(A);
  eigenvecs_exact = tridiag_sym_toeplitz_exact_eigenvecs(A);
  sort_eigen(eigenvals_exact, eigenvecs_exact);

  eigenvals_numerical = jacobi_solver(A, S, 1.E-9);
  sort_eigen(eigenvals_numerical, S);
  U1 = S * eigenvecs_exact;
  
  U2 = eigenvecs_exact;
  jacobi_solver(A, U2, 1.E-9);
  
  for (int i = 0; i < n; i++)
    REQUIRE(eigenvals_numerical[i] == Approx(eigenvals_exact[i]).epsilon(TOL));
  
  for (int i = 0; i < n; i++) {
    colvec u1 = eigenvecs_exact.col(i);
    colvec n = eigenvecs_numerical.col(i);
    REQUIRE(fabs(dot(e,n)) == Approx(norm(e)*norm(n)).epsilon(TOL));
  }
}

