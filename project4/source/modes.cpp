#include "utils.hpp"

/* Loops over temperature range T1 to T2 with steps dT. Lattice size LxL. 
   Returns result arma::mat. Simulation runs for equilibration_cycles before 
   logging of results begins. Total runs = equilibration_cycles + cycles. */
mat temp_range(int L, double entropy, double T1, double dT, double T2,
	       int cycles, int equilibration_cycles) {
  mat results;
  int n, num_accepted;

  // Setting up results matrix
  n = (T2 - T1)/dT + 1;
  if (n > MAX_SAMPLES) {
    cerr << "error: too many steps dT" << endl;
    exit(2);
  }
  results = mat(n, 6); // T, E, CV, Mabs, Chi, M

  #pragma omp parallel
  {
    double E, M, T;
    double Eacc, E2acc, Macc, M2acc, Macc_abs;
    char *p, **lattice;
    double wij[17];
    mt19937_64 rng(omp_get_thread_num());
    uniform_real_distribution<double> dist(0.,1.);

    // Set up lattice structure
    p = (char *) calloc(L*L, sizeof(char));
    CHECK_ALLOC(p);
    lattice = (char **) calloc(L, sizeof(char *));
    CHECK_ALLOC(lattice);
    for (int i = 0; i < L; i++)
      lattice[i] = &p[i*L];

    // Loop over temperature range
    #pragma omp for
    for (int i = 0; i < n; i++) {
      // Initializing
      T = T1 + i*dT;
      init_metropolis(wij, T);
      init_spins(L, entropy, lattice);

      // Monte Carlo simulation
      for (int c = 0; c < equilibration_cycles; c++)  // thermalizing...
	metropolis(L, lattice, wij, E, M, num_accepted, dist, rng);
	
      E = energy(L, lattice);
      M = magnetic_moment(L, lattice);
      Eacc = E2acc = Macc = M2acc = Macc_abs = 0.;
      for (int c = 0; c < cycles; c++) {
	Eacc += E;
	E2acc += E*E;
	Macc += M;
	M2acc += M*M;
	Macc_abs += fabs(M);
	metropolis(L, lattice, wij, E, M, num_accepted, dist, rng);
      }
      results(i,0) = T;
      results(i,1) = sample_mean(Eacc, cycles+1, L);
      results(i,2) = sample_var(E2acc, Eacc, cycles+1, L)/(T*T);
      results(i,3) = sample_mean(Macc_abs, cycles+1, L);
      results(i,4) = sample_var(M2acc, Macc_abs, cycles+1, L)/T;
      results(i,5) = sample_mean(Macc, cycles+1, L);
    }
      
    // Clean up
    free(p);
    free(lattice);
  }

  return results;
}

/* Runs Monte Carlo simulations over LxL lattice at temp. T for given number of cycles. 
   Returns result arma::mat */
mat equilibration(int L, double entropy, double T, int cycles) {
  mat results;
  int n, mod, num_accepted;
  double E, M;
  double Eacc, E2acc, Macc, M2acc, Macc_abs;
  char *p, **lattice;
  double wij[17];
  mt19937_64 rng(0);
  uniform_real_distribution<double> dist(0.,1.);

  // Setting up results matrix
  if (cycles < MAX_SAMPLES) {
    n = cycles;
    mod = 1;
  }
  else {
    n = MAX_SAMPLES;
    mod = cycles/MAX_SAMPLES;
    if (cycles % MAX_SAMPLES)
      mod++;
  }
  results = mat(n, 7); // c, E, CV, Mabs, Chi, M, acc

  // Set up lattice structure
  p = (char *) calloc(L*L, sizeof(char));
  CHECK_ALLOC(p);
  lattice = (char **) calloc(L, sizeof(char *));
  CHECK_ALLOC(lattice);
  for (int i = 0; i < L; i++)
    lattice[i] = &p[i*L];

  // Initializing
  init_metropolis(wij, T);
  init_spins(L, entropy, lattice);

  // Monte Carlo simulation
  E = energy(L, lattice);
  M = magnetic_moment(L, lattice);
  Eacc = E2acc = Macc = M2acc = Macc_abs = 0.;
  int k = 0;
  for (int c = 0; c < cycles; c++) {
    Eacc += E;
    E2acc += E*E;
    Macc += M;
    M2acc += M*M;
    Macc_abs += fabs(M);
    if (c % mod == 0) {
      results(k,0) = c;
      results(k,1) = sample_mean(Eacc, c+1, L);
      results(k,2) = sample_var(E2acc, Eacc, c+1, L)/(T*T);
      results(k,3) = sample_mean(Macc_abs, c+1, L);
      results(k,4) = sample_var(M2acc, Macc_abs, c+1, L)/T;
      results(k,5) = sample_mean(Macc, c+1, L);
      results(k,6) = num_accepted;
      k++;
    }
    num_accepted = 0;
    metropolis(L, lattice, wij, E, M, num_accepted, dist, rng);
  }
  
  // Clean up
  free(p);
  free(lattice);
  
  return results;
}


/* Estimation of the probability distribution at temperature T. 
   Returns result arma::mat with bins. Simulation runs for equilibration_cycles before 
   logging of results begins. Total runs = equilibration_cycles + cycles. */
mat probability_distribution(int L, double entropy, double T,
			     int cycles, int equilibration_cycles) {
  mat results;
  int num_accepted, idx;
  double E, M;
  double Emin = -2., Emax = 0.;
  double dE;
  char *p, **lattice;
  double wij[17];
  mt19937_64 rng(0);
  uniform_real_distribution<double> dist(0.,1.);

  // Setting up results matrix
  results = mat(PROB_DIST_BINS, 2); // bins in the normalized energy range [-2,2]
  results.zeros();
  dE = (Emax - Emin)/PROB_DIST_BINS;  // step in the normalized energy range [-2,2]
  for (idx = 0; idx < PROB_DIST_BINS; idx++)
    results(idx, 0) = Emin + idx*dE;
  
  // Set up lattice structure
  p = (char *) calloc(L*L, sizeof(char));
  CHECK_ALLOC(p);
  lattice = (char **) calloc(L, sizeof(char *));
  CHECK_ALLOC(lattice);
  for (int i = 0; i < L; i++)
    lattice[i] = &p[i*L];

  // Initializing  
  init_metropolis(wij, T);
  init_spins(L, entropy, lattice);

  // Monte Carlo simulation
  for (int c = 0; c < equilibration_cycles; c++)
    metropolis(L, lattice, wij, E, M, num_accepted, dist, rng);

  E = energy(L, lattice);
  for (int c = 0; c < cycles; c++) {
    idx = (int) ((E/(L*L) - Emin)/dE);
    if (idx >= PROB_DIST_BINS)
      idx = PROB_DIST_BINS - 1;
    results(idx,1) += 1;
    metropolis(L, lattice, wij, E, M, num_accepted, dist, rng);
  }
  
  // Clean up
  free(p);
  free(lattice);
  
  return results;
}
