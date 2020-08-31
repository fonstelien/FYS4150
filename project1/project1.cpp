#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;

colvec solve_general(rowvec a, rowvec b, rowvec c, colvec b_tilde);
colvec solve_2nd_derivative(colvec b_tilde);
colvec solve_arma(rowvec a, rowvec b, rowvec c, colvec b_tilde);


int main(int argc, char *argv[])
{
  double start_rng = 0.;
  double end_rng = 1.;
  int n = (int) 1E1;
  double h = (end_rng-start_rng)/(n+1);

  // Init diagonals
  rowvec a = randu<rowvec>(n);
  rowvec b = randu<rowvec>(n);
  rowvec c = randu<rowvec>(n);

  // Init x,b in Ax = b
  colvec x = linspace<colvec>(start_rng+h, end_rng-h, n);  // Skip range endpoints
  colvec b_tilde = 100*exp(-10*x);

  colvec solution_arma = solve_arma(a, b, c, b_tilde);
  colvec solution_general = solve_general(a, b, c, b_tilde);
  (solution_general-solution_arma).print("control general");

  rowvec base = zeros<rowvec>(n);
  solution_arma = solve_arma(base-1, base+2, base-1, b_tilde);
  colvec solution_2n_derivative = solve_2nd_derivative(b_tilde);  
  (solution_2n_derivative-solution_arma).print("control 2nd derivative");
  
  return 0;
}

colvec solve_general(rowvec a, rowvec b, rowvec c, colvec b_tilde) {
  int n = b_tilde.n_rows;
  colvec solution = b_tilde;
  rowvec b_mod = b;
  
  // Forward substitution
  for (int i = 1; i < n; i++) {
    double fac = a[i]/b_mod[i-1];
    b_mod[i] -= fac*c[i-1];
    solution[i] -= fac*solution[i-1];
  }

  // Solving for the last row
  solution[n-1] /= b_mod[n-1];

  // Backward substitution gives solution
  for (int i = n-2; i > -1; i--)
    solution[i] = (solution[i] - c[i]*solution[i+1])/b_mod[i];
  
  return solution;
}

colvec solve_2nd_derivative(colvec b_tilde) {
  int n = b_tilde.n_rows;
  colvec solution = b_tilde;

  // Forward substitution
  for (int i = 1; i < n; i++)
    solution[i] = (i+1)*solution[i] + solution[i-1];

  // Solving for the last row
  solution[n-1] /= (n+1);

  // Backward substitution gives solution
  for (int i = n-2; i > -1; i--)
    solution[i] = (solution[i] + (i+1)*solution[i+1])/(i+2);

  return solution;
}


colvec solve_arma(rowvec a, rowvec b, rowvec c, colvec b_tilde) {
  int n = b_tilde.n_rows;
  rowvec aa = a;
  rowvec cc = c;
  mat A(n,n);

  A.zeros();
  aa.shed_col(0);
  A.diag(-1) = aa;
  A.diag(0) = b;
  cc.shed_col(n-1);
  A.diag(1) = cc;

  return solve(A, b_tilde);
}
