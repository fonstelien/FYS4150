#include "utils.hpp"

/* Loops over temperature range T1 to T2 with steps dT. Lattice size LxL. 
   Returns result arma::mat. Simulation runs for equilibration_cycles before 
   logging of results begins. Total runs = equilibration_cycles + cycles. */
mat temp_range(int L, double entropy, double T1, double dT, double T2,
	       int cycles, int equilibration_cycles) {
  mat results;
  int n, num_accepted;

  // Setting up results matrix
  n = (int) ((T2 - T1)/dT + 1.);
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
      results(i,2) = sample_var(E2acc, Eacc, cycles+1, L)/(T*T);  // heat capacity CV
      results(i,3) = sample_mean(Macc_abs, cycles+1, L);
      results(i,4) = sample_var(M2acc, Macc_abs, cycles+1, L)/T;  // susceptibiliy Chi
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
mat equilibration(int L, double entropy, double T, int cycles, int equilibration_cycles) {
  mat results;
  int n, k, mod, num_accepted;
  double E, M;
  double Eacc, E2acc, Macc, M2acc, Macc_abs;
  char *p, **lattice;
  double wij[17];
  mt19937_64 rng(RNG_SEED);
  uniform_real_distribution<double> dist(0.,1.);

  // Setting up results matrix
  if (cycles+equilibration_cycles < MAX_SAMPLES) {
    n = cycles + equilibration_cycles;
    mod = 1;
  }
  else {
    n = MAX_SAMPLES;
    mod = (cycles+equilibration_cycles)/MAX_SAMPLES;
    if ((cycles+equilibration_cycles) % MAX_SAMPLES)
      mod++;
  }
  results = mat(n, 8); // c, E, CV, Mabs, Chi, M, acc, Eraw
  results.zeros();

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
  E = energy(L, lattice);
  M = magnetic_moment(L, lattice);

  // Monte Carlo simulation
  k = 0;
  for (int c = 0; c < equilibration_cycles; c++) {  // thermalizing...
    if (c % mod == 0) {
      results(k,0) = c;
      results(k,6) = num_accepted;
      results(k,7) = E/(L*L);
      k++;
    }
    metropolis(L, lattice, wij, E, M, num_accepted, dist, rng);
  }

  Eacc = E2acc = Macc = M2acc = Macc_abs = 0.;
  for (int c = 0; c < cycles; c++) {
    Eacc += E;
    E2acc += E*E;
    Macc += M;
    M2acc += M*M;
    Macc_abs += fabs(M);
    if ((c+equilibration_cycles) % mod == 0) {
      results(k,0) = c + equilibration_cycles;
      results(k,1) = sample_mean(Eacc, c+1, L);
      results(k,2) = sample_var(E2acc, Eacc, c+1, L)/(T*T);
      results(k,3) = sample_mean(Macc_abs, c+1, L);
      results(k,4) = sample_var(M2acc, Macc_abs, c+1, L)/T;
      results(k,5) = sample_mean(Macc, c+1, L);
      results(k,6) = num_accepted;
      results(k,7) = E/(L*L);
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
			     int cycles, int equilibration_cycles, int bins, double &Evar) {
  mat results;
  int num_accepted, idx;
  double E, M;
  double Eacc, E2acc;
  double Emin = -2., Emax = 2.;
  double dE;
  char *p, **lattice;
  double wij[17];
  mt19937_64 rng(RNG_SEED);
  uniform_real_distribution<double> dist(0.,1.);

  // Setting up results matrix
  results = mat(bins, 3); // bins in the normalized energy range [-2,2]
  results.zeros();
  dE = (Emax - Emin)/bins;  // step in the normalized energy range [-2,2]
  for (idx = 0; idx < bins; idx++)
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
    metropolis(L, lattice, wij, E, M, num_accepted, dist, rng);  // thermalizing...

  E = energy(L, lattice);
  Eacc = E2acc = 0.;
  for (int c = 0; c < cycles; c++) {
    Eacc += E;
    E2acc += E*E;
    idx = (int) ((E/(L*L) - Emin)/dE);
    if (idx >= bins)
      idx = bins - 1;
    results(idx,1) += 1;
    results(idx,2) += 1./cycles;
    metropolis(L, lattice, wij, E, M, num_accepted, dist, rng);
  }
  
  // Clean up
  free(p);
  free(lattice);
  
  Evar = sample_var(E2acc, Eacc, cycles+1, L);  
  return results;
}


/* Returns result arma::mat with lattice at end of simulation. */
mat get_lattice(int L, double entropy, double T, int cycles) {
  mat results;
  int num_accepted;
  double E, M;
  char *p, **lattice;
  double wij[17];
  mt19937_64 rng(RNG_SEED);
  uniform_real_distribution<double> dist(0.,1.);

  // Setting up results matrix
  results = mat(L, L);
  
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
  for (int c = 0; c < cycles; c++)
    metropolis(L, lattice, wij, E, M, num_accepted, dist, rng);

  for (int i = 0; i < L; i++)
    for (int j = 0; j < L; j++)
      results(i,j) = lattice[i][j];
  
  // Clean up
  free(p);
  free(lattice);
  
  return results;
}


