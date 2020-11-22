#include <iostream>
#include <cmath>
#include <random>
#include <armadillo>
#include <unistd.h>
#include <omp.h>

// RNG seed
#define RNG_SEED 0

// Maximum number of stored samples in the results matrix
#define MAX_SAMPLES 1000

#define DEBUG(msg) cout << "DEBUG " << msg << endl;

// Checks if memory allocation was successful.
#define CHECK_ALLOC(allocation) {					\
    if(allocation == NULL) {						\
      fprintf(stderr, "error: could not allocate memory.\n");		\
      exit(5);}}

// The beta factor in the Boltzmann distribution. Scale with kB=1.
#define BETA(T) (1/T)

// Periodic boundary indexes
// NOTE: L is the lattice side length and must be declared exactly as L in the namespace
#define PERIOD(idx) ((idx + L) % L)

using namespace std;
using namespace arma;

/* From utils.cpp */

/* Initializes lattice. Spins are either 1 (up) or -1 (down). Entropy s/2 is the 
   likelihood of any spin -1 such that s=0.0 gives all 1s; s=1.0 gives complete 
   random spins. */
void init_spins(int L, double s, char **lattice);

/* Returns the energy/J of the lattice. */
double energy(int L, char **lattice);

/* Returns the magnetic moment of the lattice. */
double magnetic_moment(int L, char **lattice);

/* Initializes metropolis acceptance rule vector with the Boltzmann distribution 
   for the energy change from state i->j at temp. T */
void init_metropolis(double *wij, double T);

/* Runs the metropolis algorithm with monte carlo selection on the lattice. Updates
   total energy E/J and total magnetic moment along the way. */
void metropolis(int L, char **lattice, double *wij, double &E, double &M, int &num_accepted,
		uniform_real_distribution<double> &dist, mt19937_64 &rng);

/* Calculates the n-sample sample-mean of variable x, normalized to lattice size N=L*L. */
double sample_mean(double x, int n, int L);

/* Calculates the n-sample sample-var of variable x, normalized to lattice size N=L*L. */
double sample_var(double x2, double x, int n, int L);

/* Plots the lattice to stdout. */
void plot_lattice(int L, char **lattice);


/* From modes.cpp */

/* Loops over temperature range T1 to T2 with steps dT. Lattice size LxL. 
   Returns result arma::mat. Simulation runs for equilibration_cycles before 
   logging of results begins. Total runs = equilibration_cycles + cycles. */
mat temp_range(int L, double entropy, double T1, double dT, double T2,
	       int cycles, int equilibration_cycles);

/* Runs Monte Carlo simulations over LxL lattice at temp. T for given number of cycles. 
   Returns result arma::mat */
mat equilibration(int L, double entropy, double T, int cycles);

/* Estimation of the probability distribution at temperature T. 
   Returns result arma::mat with bins. Simulation runs for equilibration_cycles before 
   logging of results begins. Total runs = equilibration_cycles + cycles. */
mat probability_distribution(int L, double entropy, double T,
			     int cycles, int equilibration_cycles, int bins, double &Evar);

/* Returns result arma::mat with lattice at end of simulation. */
mat get_lattice(int L, double entropy, double T, int cycles);
