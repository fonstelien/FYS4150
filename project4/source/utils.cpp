#include "utils.hpp"

using namespace std;

/* Initializes lattice with periodic boundaries. Spins are either 1 (up) or 
   -1 (down). Entropy s/2 is the likelihood of any spin -1 such that s=0.0 gives all 1s;
   s=1.0 gives complete random spins. */
void init_spins(int L, double s, char **lattice) {
  for (int i = 0; i < L; i++)
    for (int j = 0; j < L; j++) {
      lattice[i][j] = 1;
      if ((double) rand()/RAND_MAX < s*.5)
	lattice[i][j] = -1;
    }
}

/* Returns the total energy/J of the lattice. */
double energy(int L, char **lattice) {
  int e = 0;

  for (int i = 0; i < L; i++)
    for (int j = 0; j < L; j++)
      e -= (int) lattice[i][j]*(lattice[i][PERIOD(j+1)] + \
				lattice[PERIOD(i+1)][j]);

  return (double) e;
}

/* Returns the total magnetic moment of the lattice. */
double magnetic_moment(int L, char **lattice) {
  int m = 0;

  for (int i = 0; i < L; i++)
    for (int j = 0; j < L; j++)
      m += (int) lattice[i][j];

  return (double) m;
}

/* Initializes metropolis acceptance rule vector with the Boltzmann distribution 
   for the energy change from state i->j at temp. T */
void init_metropolis(double *wij, double T) {
  wij[0] = std::exp(-BETA(T)*(-8));
  wij[4] = std::exp(-BETA(T)*(-4));
  wij[8] = 1.;
  wij[12] = std::exp(-BETA(T)*4);
  wij[16] = std::exp(-BETA(T)*8);
}

/* Runs the metropolis algorithm with monte carlo selection on the lattice. Updates
   total energy E/J and total magnetic moment along the way. */
void metropolis(int L, char **lattice, double *wij, double &E, double &M,
		uniform_real_distribution<double> &dist, mt19937_64 &rng) {
  int x,y;
  int dE;
  int b;
  
  for (int i = 0; i < L; i++)
    for (int j = 0; j  < L; j++) {
      x = dist(rng)*L;
      y = dist(rng)*L;
      dE = 2*lattice[x][y]*(lattice[PERIOD(x-1)][y] +	\
			    lattice[PERIOD(x+1)][y] +	\
			    lattice[x][PERIOD(y-1)] +   \
			    lattice[x][PERIOD(y+1)]);
      b = (int) (dist(rng) < wij[dE+8]);
      lattice[x][y] -= 2*b*lattice[x][y];
      E += b*dE;
      M += 2*b*lattice[x][y];
    }
}

/* Calculates the n-sample sample-mean of variable x, normalized to lattice size N=L*L. */
double sample_mean(double x, int n, int L) {
  int N = L*L;
  return x/n/N;
}

/* Calculates the n-sample sample-var of variable x, normalized to lattice size N=L*L. */
double sample_var(double x2, double x, int n, int L) {
  int N = L*L;
  double mux = x/n;
  return (x2/n - mux*mux)/N;
}


/* Plots the lattice to stdout. */
void plot_lattice(int L, char **lattice) {
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      if (lattice[i][j] > 0)
	cout << "+ ";
      else
	cout << "  ";
    }
    cout << endl;
  }
}
