#include "utils.hpp"

using namespace std;

// Modes
#define TEMP_RANGE 1
#define EQUILIBRATION 2
#define PROBABILITY 3

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

  // Parsing args
  if (argc < 2) {
    cerr << "error: missing arguments" << endl;
    print_usage();
    exit(1);
  }

  opterr = 0;
  while ((opt = getopt(argc, argv, "hRE:P:pl:t:c:e:s:")) != -1) {
    switch (opt) {
    case 'h':
      print_usage();
      exit(0);
    case 'R':
      mode = TEMP_RANGE;
      T1 = strtod(optarg, NULL);
      continue;
    case 'E':
      mode = EQUILIBRATION;
      T1 = strtod(optarg, NULL);
      continue;
    case 'P':
      mode = PROBABILITY;
      T1 = strtod(optarg, NULL);
      continue;
    // case 'p':
    //   print_lattice = true;
    //   continue;
    case 'l':
      L = atoi(optarg);
      continue;
    case 't':
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
    results = probability_distribution(L, entropy, T1, cycles, equilibration_cycles);
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
    cout << "n,E" << endl;
    results.save(cout, csv_ascii);
    break;
    
  default:
    break;
  }

  return 0;
}
