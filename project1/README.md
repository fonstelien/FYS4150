# FYS4150 H20 Project 1
You will find my implementation of the algorithms in project1.cpp. Compile with
```
$ make
g++ -std=c++11 project1.cpp -o project1 -Wall -Wextra -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result
```
... and then have a look at the `usage` for instructions on how to run it. Enjoy! :whale:
```
$ ./project1 -h
usage: project1 [-h | -a | -g | -d | -l | -t] [-epc] n
  -h   print this text
  -a   run all
  -g   run general tridiagonal
  -d   run optimized for 2nd derivative
  -l   run using armadillo LU decomposition
  -t   test all implementations against arma::solve()
       exit value 0 indicates success;
       11,12,13 failure in general, 2nd derivative, or LU decomposition
  -e   print log10 of relative error
  -p   print numerical solution to stdout in csv format (see Examples)
  -c   print closed-form solution to stdout in csv format (see Examples)
   n   number of calculation points

 Results write to stdout: 1st pos is always n.
 For -a option CPU time are in order general, 2nd deriv., LU decomp.

 Examples:
 $ project1 -a 1000
 1000, 0.492, 0.131, 422.61
 $ project1 -dp 4
 xi, v[xi]
 2.00000e-01, 1.44596e-01
 4.00000e-01, 2.87850e-01
 6.00000e-01, 4.21188e-01
 8.00000e-01, 4.81265e-01

```
