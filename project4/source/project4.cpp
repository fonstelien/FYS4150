#include "utils.hpp"

using namespace std;

// Modes
#define TEMP_RANGE 1
#define EQUILIBRATION 2

// Maximum number of stored samples
#define MAX_SAMPLES 1000

/* Prints program usage */
void print_usage() {
  cout << "usage: project3 [-h | -f | -E | -V | -T | -C | -P] -n -y -b -v -j -s -r" << endl;
  cout << "..." << endl;
}

int main(int argc, char *argv[]) {
  // options
  int opt;
  int mode = TEMP_RANGE;
  int L = 20;
  double T1 = 2.0, dT = .1, T2 = 2.4;
  double entropy = 1.0;
  int cycles = (int) 1.E5;
  int equilibration_cycles = 0;
  int n = -1;
  int mod = -1;
  
  // other vars
  char *p;
  char **lattice = &p;
  double wij[17];
  double E, M, T;
  double Eacc, E2acc, Macc, M2acc, Macc_abs;
  double energies[MAX_SAMPLES], magnetic_moments[MAX_SAMPLES], abs_moments[MAX_SAMPLES];
  double x[MAX_SAMPLES];
  mat results;
  bool print_lattice = false;

  // Parsing args
  if (argc < 2) {
    cerr << "error: missing arguments" << endl;
    print_usage();
    exit(1);
  }

  opterr = 0;
  while ((opt = getopt(argc, argv, "hpE:L:T:c:e:s:")) != -1) {
    switch (opt) {
    case 'h':
      print_usage();
      exit(0);
    case 'p':
      print_lattice = true;
      continue;
    case 'E':
      mode = EQUILIBRATION;
      T1 = strtod(optarg, NULL);
      continue;
    case 'L':
      L = atoi(optarg);
      continue;
    case 'T':
      T1 = strtod(optarg, NULL);
      dT = strtod(argv[optind++], NULL);
      T2 = strtod(argv[optind], NULL);
      continue;
    case 'c':
      cycles = (int) strtod(optarg, NULL);
      continue;
    case 'e':
      equilibration_cycles = (int) strtod(optarg, NULL);
      continue;
    case 's':
      entropy = strtod(optarg, NULL);
      continue;
    default:
      cerr << "error: unknown option: " << (char) optopt << endl;
      print_usage();
      exit(1);
    }
  }

  // Post-process parsed args
  switch (mode) {
  case TEMP_RANGE:
    n = (T2 - T1)/dT + 1;
    if (n > MAX_SAMPLES) {
      cerr << "error: too many steps dT" << endl;
      exit(2);
    }
    results = mat(n, 6); // T, E, CV, Mabs, Chi, M
    break;

  case EQUILIBRATION:
    if (cycles < MAX_SAMPLES)
      n = cycles;
    else {
      n = MAX_SAMPLES;
      mod = cycles/MAX_SAMPLES;
      if (cycles % MAX_SAMPLES)
	mod++;
    }
    results = mat(n, 6); // c, E, CV, Mabs, Chi, M
    break;
  }
  
  // Run simulations
  switch (mode) {
  case TEMP_RANGE:
    #pragma omp parallel private(p, lattice, wij, E, M, T)
    {
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
	T = T1 + i*dT;
	init_metropolis(wij, T);
	init_spins(L, entropy, lattice);

	// Monte Carlo simulation
	for (int c = 0; c < equilibration_cycles; c++)  // thermalizing...
	  metropolis(L, lattice, wij, E, M, dist, rng);
	
	E = energy(L, lattice);
	M = magnetic_moment(L, lattice);
	Eacc = E2acc = Macc = M2acc = Macc_abs = 0.;
	for (int c = 0; c < cycles; c++) {
	  Eacc += E;
	  E2acc += E*E;
	  Macc += M;
	  M2acc += M*M;
	  Macc_abs += fabs(M);
	  metropolis(L, lattice, wij, E, M, dist, rng);
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
    break;

  case EQUILIBRATION:
    {
      mt19937_64 rng(0);
      uniform_real_distribution<double> dist(0.,1.);

      // Set up lattice structure
      p = (char *) calloc(L*L, sizeof(char));
      CHECK_ALLOC(p);
      lattice = (char **) calloc(L, sizeof(char *));
      CHECK_ALLOC(lattice);
      for (int i = 0; i < L; i++)
	lattice[i] = &p[i*L];

      init_metropolis(wij, T1);
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
	  results(k,2) = sample_var(E2acc, Eacc, c+1, L)/(T1*T1);
	  results(k,3) = sample_mean(Macc_abs, c+1, L);
	  results(k,4) = sample_var(M2acc, Macc_abs, c+1, L)/T1;
	  results(k,5) = sample_mean(Macc, c+1, L);
	  k++;
	}
	metropolis(L, lattice, wij, E, M, dist, rng);
      }
      
    }
    break;

  default:
    break;
  }


  // Print output
  cout << "# $";
  for (int i = 0; i < argc; i++)
    cout << " " << argv[i];
  cout << endl;

  switch (mode) {
  case TEMP_RANGE:
      cout << "T,E,CV,Mabs,Chi,M" << endl;
      results.save(cout, csv_ascii);
    break;

  case EQUILIBRATION:
    if (print_lattice)
      plot_lattice(L, lattice);
    else {
      cout << "c,E,CV,Mabs,Chi,M" << endl;
      results.save(cout, csv_ascii);
    }
    break;

  default:
    break;
  }

  // Clean up
  if (p)
    free(p);
  if (lattice)
    free(lattice);
  
  return 0;
}
