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
  cout << "usage: project4 [-h | R | E | P | L ] -lcesb" << endl;
  cout << " -h  prints this text" << endl;
  cout << " -R  runs over temperature Range: '-R T1 dT T2'" << endl;
  cout << " -E  runs Equilibration at given temp: -E T" << endl;
  cout << " -P  Probability distribution estimation at given temp: -P T" << endl;
  cout << " -L  outputs the Lattice at end of sim. for given temp: -L T" << endl;
  cout << " -l  side length of lxl lattice. Default l=20" << endl;
  cout << " -c  cycles in simulation" << endl;
  cout << " -e  equilibration cycles for opts -R and -P. Default e=c/10" << endl;
  cout << " -s  initial entropy in range [0.0,1.0]" << endl;
  cout << " -b  number of bins in the prob. density estimation for -P. Default b=20" << endl;
  cout << endl;
  cout << " All output in csv format." << endl;
  cout << " Example:" << endl;
  cout << " $ project4 -R 2 .1 2.3 -l 20 -c 1e5" << endl;
  cout << " T,E,CV,Mabs,Chi,M" << endl;
  cout << " 2.0000000000000000e+00,-1.7459954400455997e+00,7.3134602046044161e-01,9.1147373526264741e-01,3.8114415508978711e-01,-9.1147373526264741e-01" << endl; 
  cout << " 2.1000000000000001e+00,-1.6617241827581726e+00,9.6394348202441926e-01,8.6970120298797016e-01,8.1059688821260345e-01,-8.6970120298797016e-01" << endl;
  cout << " 2.2000000000000002e+00,-1.5496672033279668e+00,1.3470200211441872e+00,7.8916300836991626e-01,3.2399410438633409e+00,-1.3165198348016521e-01" << endl;
}

int main(int argc, char *argv[]) {
  int opt;
  int mode = -1;
  int L = 20;
  double T1 = 2.0, dT = .1, T2 = 2.4;
  double entropy = 1.0;
  int cycles = (int) 1.E5;
  int equilibration_cycles = -1;
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
  while ((opt = getopt(argc, argv, "hR:E:P:L:l:c:e:s:b:")) != -1) {
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

  if (equilibration_cycles < 0)
    equilibration_cycles = cycles/10;

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
