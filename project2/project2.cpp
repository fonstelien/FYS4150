#include <fstream>
#include <unistd.h>
#include <time.h>
#include "utils.hpp"

#define TOL 1.E-3

using namespace arma;
using namespace std;

/* From utils.cpp */


int main(int argc, const char **argv) {
  vec eigenvals_numerical, eigenvals_exact;
  mat A, S;
  double rho_max, eps;
  int n;
  float multiplier = 0;
  ofstream out_file;
  clock_t start_t, end_t;
  double lapsed_time;
  
  /* Parsing args */
  rho_max = strtod(argv[1], NULL);  
  n = atoi(argv[2]);
  multiplier = strtod(argv[3], NULL);
  out_file.open(argv[4]);

  
  eigenvals_exact = vec(n);
  eigenvals_exact[0] = 3.;
  for (int i = 1; i < n; i++)
    eigenvals_exact[i] = eigenvals_exact[i-1] + 4.;

  eps = TOL + 1;
  for (; eps > TOL; n *= multiplier) {
    A = tridiag_sym_general(n, rho_max);
    S.eye(n,n);

    start_t = clock();
    eigenvals_numerical = jacobi_solver(A, S, 1.E-3);
    end_t = clock();
    
    sort_eigen(eigenvals_numerical, S);
    eps = fabs((eigenvals_exact[0] - eigenvals_numerical[0])/eigenvals_exact[0]);

    lapsed_time = (double) (end_t - start_t) / CLOCKS_PER_SEC;
    cout << "n, secs, eps: " << n << ", " << lapsed_time << ", " << eps << endl;
  }
  
  for (int i = 0; i < (3 < n ? 3 : n); i++)
    cout << eigenvals_numerical[i] << " " << eigenvals_exact[i] << endl;

  out_file << n << endl;
  out_file << eigenvals_numerical;
  out_file << S << endl;
  out_file.close();

  
  return 0;
}
