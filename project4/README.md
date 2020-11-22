# FYS4150 H20 Project 4
You will find the report and supporting files are in the `report` folder and source code in the `source` folder. The source code is arranged such that you select which mode to run the main program in `project4.cpp` in, and these modes are implemented in the `modes.cpp` file. Other stuff that is common to the modes are implemented in `utils.cpp`, and `utils.hpp` ties it all together.

Parallelization of the `temp_range()` function is done in OpenMP, which is my preferred parallelization environment, and which I also think is the right choice since I am running all simulations on my laptop.

All `*.csv`s used for making the plots are found in the `results` folder, together with my `project4.ipynb`. Most of the results are easily reproducible with the accompanying `bash` scripts. You can also play with the Jupyter Notebook file and I have attached a `fys4150.yml` which you can use to clone my environment to your own machine.

To compile the C++ code, move to `source` and do
```
$ make
g++ -std=c++11 -c utils.cpp -o utils.o -Wall -Wextra -fopenmp -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result 
g++ -std=c++11 -c modes.cpp -o modes.o -Wall -Wextra -fopenmp -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result 
g++ -std=c++11 project4.cpp utils.o modes.o -o project4 -Wall -Wextra -fopenmp -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result
```
... and then have a look at the `usage` for instructions on how to run it. Enjoy! :whale:
```
$ ./project4 -h
usage: project4 [-h | R | E | P | L ] -lcesb
 -h  prints this text
 -R  runs over temperature Range: '-R T1 dt T2'
 -E  runs Equilibration at given temp: -E T
 -P  Probability distribution estimation at given temp: -P T
 -L  outputs the Lattice at end of sim. for given temp: -L T
 -l  side length of lxl lattice. Default l=20
 -c  cycles in simulation
 -e  equilibration cycles for opts -R and -P. Default e=c/10
 -s  initial entropy in range [0.0,1.0]
 -b  number of bins in the prob. density estimation for -P. Default b=20

 All output in csv format.
 Example:
 $ project4 -R 2 .1 2.3 -l 20 -c 1e5
 T,E,CV,Mabs,Chi,M
 2.0000000000000000e+00,-1.7459954400455997e+00,7.3134602046044161e-01,9.1147373526264741e-01,3.8114415508978711e-01,-9.1147373526264741e-01
 2.1000000000000001e+00,-1.6617241827581726e+00,9.6394348202441926e-01,8.6970120298797016e-01,8.1059688821260345e-01,-8.6970120298797016e-01
 2.2000000000000002e+00,-1.5496672033279668e+00,1.3470200211441872e+00,7.8916300836991626e-01,3.2399410438633409e+00,-1.3165198348016521e-01
```