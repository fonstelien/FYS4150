# FYS4150 H20 Project 2
You will find all source code in the `soruce` folder. Compile with
```
source$ make
g++ -std=c++11 -c utils.cpp -o utils.o -Wall -Wextra -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result 
g++ -std=c++11 -c jacobi.cpp -o jacobi.o -Wall -Wextra -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result 
g++ -std=c++11 -c polynomial_expansion.cpp -o polynomial_expansion.o -Wall -Wextra -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result 
g++ -std=c++11 project2.cpp utils.o jacobi.o polynomial_expansion.o -o project2 -Wall -Wextra -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result
g++ -std=c++11 -c test_define.cpp -o test_define.o -Wall -Wextra -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result 
g++ -std=c++11 test_project2.cpp utils.o jacobi.o polynomial_expansion.o test_define.o -o test_project2 -Wall -Wextra -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result
```
... and then have a look at the `usage` for instructions on how to run it. Unit tests are found in `test_project2.cpp`, and you run it with `source$ ./test_project2` Enjoy! :whale:
```
source$ ./project2 -h
usage: project2 [-h | -C | -I | -V | -T] [-j | -p] [-t | -g] [-re] n
  -h   print this text
  -C   Compare exact and numerical results. Prints to stdout.
  -I   number of Iterations before convergence (Jacobi only)
  -V   prints eigenVectors to stdout (Jacobi only)
  -T   prints CPU Time for algorithm
  -j   run Jacobi algorithm
  -p   run Polynomial expansion algorithm
  -t   Toeplitz tridiagonal matrix
  -g   'general' tridiagonal matrix
  -r   rho max in the 'general' matrix
  -e   convergence tolerance (eps)
  -n   
Where suitable, the results are printed in csv format.
example:
$ project2 -Cjg -e 1.E-2 -r 10.0 -n 100
exact,numeric
3.000000000000e+00,2.779213340881e+00
7.000000000000e+00,6.654696896637e+00
1.100000000000e+01,1.054835747359e+01
...
```