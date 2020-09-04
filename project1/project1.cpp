#include <iostream>
#include <armadillo>
#include <unistd.h>
#include <time.h>

using namespace arma;
using namespace std;

#define LAPSED_TIME(t0) (((double) (clock() - t0)) / CLOCKS_PER_SEC * 1000)

// Calculation range
#define START_RANGE 0.
#define END_RANGE 1.

// Program modes
#define RUN_ALL 0
#define GENERAL 1
#define DERIVATIVE 2
#define LIBRARY 3
#define EPS 4
#define CLOSED_FORM 5
#define TEST 6

/* Prints program usage */
void print_usage() {
  cout << "usage: project1 [-h | -a | -g | -d | -l | -t] [-epc] n" << endl;
  cout << "  -h   print this text" << endl;
  cout << "  -a   run all" << endl;
  cout << "  -g   run general tridiagonal" << endl;
  cout << "  -d   run optimized for 2nd derivative" << endl;
  cout << "  -l   run using armadillo LU decomposition" << endl;
  cout << "  -t   test all implementations against arma::solve()" << endl;
  cout << "       exit value 0 indicates success;" << endl;
  cout << "       11,12,13 failure in general, 2nd derivative, or LU decomposition" << endl;
  cout << "  -e   print log10 of relative error" << endl;
  cout << "  -p   print numerical solution to stdout in csv format (see Examples)" << endl;
  cout << "  -c   print closed-form solution to stdout in csv format (see Examples)" << endl;
  cout << "   n   number of calculation points" << endl;
  cout << "" << endl;
  cout << " Results write to stdout: 1st pos is always n." << endl;
  cout << " For -a option CPU time are in order general, 2nd deriv., LU decomp." << endl;
  cout << "" << endl;
  cout << " Examples:" << endl;
  cout << " $ project1 -a 1000" << endl;
  cout << " 1000, 0.492, 0.131, 422.61" << endl;
  cout << " $ project1 -dp 4" << endl;
  cout << " xi, v[xi]" << endl;
  cout << " 2.00000e-01, 1.44596e-01" << endl;
  cout << " 4.00000e-01, 2.87850e-01" << endl;
  cout << " 6.00000e-01, 4.21188e-01" << endl;
  cout << " 8.00000e-01, 4.81265e-01" << endl;
}

colvec solve_tridiag(rowvec a, rowvec b, rowvec c, colvec b_tilde);
colvec solve_2nd_derivative(colvec b_tilde);
colvec solve_LU(mat L, mat U, colvec b_tilde);
colvec closed_form_solution(colvec x);
colvec log_eps(colvec exact, colvec approximated);
mat generate_A_matrix(rowvec a, rowvec b, rowvec c);
int run_test(mat L, mat U, rowvec a, rowvec b, rowvec c, colvec b_tilde);

int main(int argc, char *argv[])
{
  int n;
  int opt;
  int mode = -1;
  bool print_solution = false;
  double h;
  double max_eps = -1;
  colvec solution;
  rowvec a, b, c;
  mat L, U;
  clock_t start_time;
  double running_times[10];
  
  // Parsing args
  if (argc < 2) {
    cerr << "error: missing arguments" << endl;
    print_usage();
    exit(1);
   }

  opterr = 0;
  while ((opt = getopt(argc, argv, "hagdltepc")) != -1) {
    switch (opt) {
    case 'h':
      print_usage();
      exit(0);
    case 'a':
      mode = RUN_ALL;
      continue;
    case 'g':
      mode = GENERAL;
      continue;
    case 'd':
      mode = DERIVATIVE;
      continue;
    case 'l':
      mode = LIBRARY;
      continue;
    case 't':
      mode = TEST;
      continue;
    case 'e':
      mode = EPS;
      continue;
    case 'p':
      print_solution = true;
      continue;
    case 'c':
      mode = CLOSED_FORM;
      print_solution = true;
      continue;
    default:
      cerr << "error: unknown option: " << (char) optopt << endl;
      print_usage();
      exit(1);
    }
  }

  if (argc - optind < 1 || mode < 0) {
    cerr << "error: missing arguments" << endl;
    print_usage();
    exit(1);
  }

  n = atoi(argv[optind++]);
  h = (END_RANGE-START_RANGE)/(n+1);
  
  // Initialize vector of constants b_tilde (b in the Ax = b equation)
  colvec x = linspace<colvec>(START_RANGE+h, END_RANGE-h, n);  // Skip range endpoints
  colvec b_tilde = flipud(h*h*100*exp(-10*x));

  if (mode != DERIVATIVE) {
    a = zeros<rowvec>(n) - 1;
    b = zeros<rowvec>(n) + 2;
    c = zeros<rowvec>(n) - 1;
  }

  if (mode == RUN_ALL || mode == LIBRARY || mode == TEST) {
    mat A;

    A = generate_A_matrix(a, b, c);
    lu(L, U, A);  // Permutation matrix not needed
  }

  
  // Run calculations
  switch (mode) {
  case RUN_ALL:
    // Gemeral
    start_time = clock();
    solution = solve_tridiag(a, b, c, b_tilde);
    running_times[GENERAL] = LAPSED_TIME(start_time);

    // 2nd derivative
    start_time = clock();
    solution = solve_2nd_derivative(b_tilde);
    running_times[DERIVATIVE] = LAPSED_TIME(start_time);

    // Library
    start_time = clock();
    solution = solve_LU(L, U, b_tilde);
    running_times[LIBRARY] = LAPSED_TIME(start_time);
    break;
    
  case GENERAL:
    start_time = clock();
    solution = solve_tridiag(a, b, c, b_tilde);
    running_times[GENERAL] = LAPSED_TIME(start_time);
    break;
    
  case DERIVATIVE:
    start_time = clock();
    solution = solve_2nd_derivative(b_tilde);
    running_times[DERIVATIVE] = LAPSED_TIME(start_time);
    break;
    
  case LIBRARY:
    start_time = clock();
    solution = solve_LU(L, U, b_tilde);
    running_times[LIBRARY] = LAPSED_TIME(start_time);
    break;

  case EPS:
    {
      solution = solve_2nd_derivative(b_tilde);
      colvec exact = closed_form_solution(x);
      colvec eps = log_eps(exact, solution);
      if (!eps.is_empty())
	max_eps = max(eps);
    }
    break;
    
  case CLOSED_FORM:
    solution = closed_form_solution(x);
    break;
    
  case TEST:
    exit(run_test(L, U, a, b, c, b_tilde));
    
  default:
    cerr << "error: unknown mode" << endl;
    exit(-1);
  }

  // Print results
  // Print calculation results
  if (print_solution) {
    cout << scientific;
    cout.precision(5);
    cout << "xi, v[xi]" << endl;
    for (int i = 0; i < n; i++)
      cout << x[i] << ", " << solution[i] << endl;
  }

  // Print CPU time and eps
  else {  
    cout.precision(5);
    cout << n << ", ";
    switch (mode) {
    case RUN_ALL:
      cout << running_times[GENERAL] << ", " << running_times[DERIVATIVE] << ", " \
	   << running_times[LIBRARY] << endl;
      break;
    case GENERAL:
      cout << running_times[GENERAL] << endl;
      break;
    case DERIVATIVE:
      cout << running_times[DERIVATIVE] << endl;
      break;
    case LIBRARY:
      cout << running_times[LIBRARY] << endl;
      break;
    case EPS:
      cout << max_eps << endl;
      break;
    }
  }
  
  return 0;
}

colvec solve_tridiag(rowvec a, rowvec b, rowvec c, colvec b_tilde) {
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

mat generate_A_matrix(rowvec a, rowvec b, rowvec c) {
  int n = a.n_cols;
  rowvec aa = a;
  rowvec cc = c;
  mat A(n,n);
  
  // Set up matrix with diagonals
  A.zeros();
  aa.shed_col(0);
  A.diag(-1) = aa;
  A.diag(0) = b;
  cc.shed_col(n-1);
  A.diag(1) = cc;

  return A;
}

colvec solve_LU(mat L, mat U, colvec b_tilde) {
  colvec y;

  y = solve(L, b_tilde);
  return solve(U, y);
}

colvec closed_form_solution(colvec x) {
  colvec solution;

  solution = flipud(1 - (1 - std::exp((double) -10.))*x - exp(-10*x));
  return solution;
}

colvec log_eps(colvec exact, colvec approximated) {
  colvec eps;

  // Check if exact solution vector has zeros
  if (min(abs(exact)) == 0.) {
    cerr << "error: cannot calculate relative error; "		\
	 << "exact solution contains one or more zeros." << endl;
    return eps;
  }

  eps = abs((approximated - exact)/exact);
  return log10(eps);
}

int run_test(mat L, mat U, rowvec a, rowvec b, rowvec c, colvec b_tilde) {
  mat A;
  colvec arma_solution;

  A = generate_A_matrix(a, b, c);
  arma_solution = solve(A, b_tilde);
  
  if (max(abs(arma_solution - solve_tridiag(a, b, c, b_tilde))) > 1E-10) {
    cerr << "test error: solve_tridiag()" << endl;
    return 11;    
  }

  if (max(abs(arma_solution - solve_2nd_derivative(b_tilde))) > 1E-10) {
    cerr << "test error: solve_2nd_derivative()" << endl;
    return 12;      
  }

  if (max(abs(arma_solution - solve_LU(L, U, b_tilde))) > 1E-10) {
    cerr << "test error: solve_LU()" << endl;
    return 13;    
  }

  return 0;
}
