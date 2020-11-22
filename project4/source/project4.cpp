#include "utils.hpp"

using namespace std;

// Modes
#define TEMP_RANGE 1
#define EQUILIBRATION 2
#define PROBABILITY 3
#define LATTICE_OUT 4

// Maximum number of stored samples
#define MAX_SAMPLES 1000

/* Prints program usage */
void print_usage() {
  cout << "usage: project3 [-h | -f | -E | -V | -T | -C | -P] -n -y -b -v -j -s -r" << endl;
  cout << "..." << endl;
}

int main(int argc, char *argv[]) {
  int opt;
  int mode = -1;
  int L = 20;
  double T1 = 2.0, dT = .1, T2 = 2.4;
  double entropy = 1.0;
  int cycles = (int) 1.E5;
  int equilibration_cycles = 0;
  mat results;
  int bins = 20;  // default number of bins for the probability distribution estimation
  double Evar = 0.;

  // Parsing args
  if (argc < 2) {
    cerr << "error: missing arguments" << endl;
    print_usage();
    exit(1);
  }

  opterr = 0;
  while ((opt = getopt(argc, argv, "hR:E:P:L:pl:c:e:s:b:")) != -1) {
    switch (opt) {
    case 'h':
      print_usage();
      exit(0);
    case 'R':
      mode = TEMP_RANGE;
      T1 = strtod(optarg, NULL);
      dT = strtod(argv[optind++], NULL);
      T2 = strtod(argv[optind], NULL);
      continue;
    case 'E':
      mode = EQUILIBRATION;
      T1 = strtod(optarg, NULL);
      continue;
    case 'P':
      mode = PROBABILITY;
      T1 = strtod(optarg, NULL);
      continue;
    case 'L':
      mode = LATTICE_OUT;
      T1 = strtod(optarg, NULL);
      continue;
    case 'l':
      L = atoi(optarg);
      continue;
    case 't':
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
    case 'b':
      bins = (int) strtod(optarg, NULL);
      continue;
    default:
      cerr << "error: unknown option: " << (char) optopt << endl;
      print_usage();
      exit(1);
    }
  }

  // Post-process parsed args  
  if (mode < 1) {
    cerr << "error: select mode"<< endl;
    print_usage();
    exit(2);
  }

  // Run simulations
  switch (mode) {
  case TEMP_RANGE:
    results = temp_range(L, entropy, T1, dT, T2, cycles, equilibration_cycles);
    break;

  case EQUILIBRATION:
    results = equilibration(L, entropy, T1, cycles);
    break;

  case PROBABILITY:
    results = probability_distribution(L, entropy, T1, cycles, equilibration_cycles,
				       bins, Evar);
    break;

  case LATTICE_OUT:
    results = get_lattice(L, entropy, T1, cycles);
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
    cout << "c,E,CV,Mabs,Chi,M,accs" << endl;
    results.save(cout, csv_ascii);
    break;

  case PROBABILITY:
    cout << "# Evar=" << Evar << endl;
    cout << "E,n,p" << endl;
    results.save(cout, csv_ascii);
    break;

  case LATTICE_OUT:
    results.save(cout, csv_ascii);
    break;
    
  default:
    break;
  }

  return 0;
}
