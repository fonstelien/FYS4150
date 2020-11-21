# FYS4150 H20 Project 3
The report and supporting files are in the `report` folder. The `.csv` files with the results used for the figures in the report were all above GitHub's 50MB limit, so you will not find them in the `results` folder. You will however will find my  `project3.ipynb` along with `fys4150.yml` which you can use to clone my environment, and I will find a way to share the raw data with you upon request.

In the `source` folder you will find all source code developed for the project. Compile with
```
source$ make
g++ -std=c++11 -c non_OO_earth_sun.cpp -o non_OO_earth_sun.o -Wall -Wextra -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result 
g++ -std=c++11 -c planet.cpp -o planet.o -Wall -Wextra -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result 
g++ -std=c++11 -c solver.cpp -o solver.o -Wall -Wextra -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result 
g++ -std=c++11 -c perihelion_of_mercury.cpp -o perihelion_of_mercury.o -Wall -Wextra -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result 
g++ -std=c++11 project3.cpp non_OO_earth_sun.o planet.o solver.o perihelion_of_mercury.o -o project3 -Wall -Wextra -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result
g++ -std=c++11 -c catch.cpp -o catch.o -Wall -Wextra -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result 
g++ -std=c++11 test_project3.cpp non_OO_earth_sun.o planet.o solver.o perihelion_of_mercury.o catch.o -o test_project3 -Wall -Wextra -larmadillo -O3 -Wno-unused -Wno-unused-parameter -Wno-unused-result
```
... and then have a look at the `usage` for instructions on how to run it. Unit tests are found in `test_project3.cpp`, and you run it with `source$ ./test_project3`. It should all run smoothly... Enjoy! :whale:
```
$ ./project3 -h
usage: project3 [-h | -f | -E | -V | -T | -C | -P] -n -y -b -v -j -s -r
  -h   Print this text
  -f   Run system from csv file input. Format: mass,x,y,z,vx,vy,vz
       mass is planet's mass relative to Sun. vx,vy,vz in AU/day.
       Outputs total kinetic,potential energy, angular momentum,
       x,y,z,vx,vy,vz,ep,ek,am, for each planet. vx,vy,vz in AU/yr.
  -E   Euler's forward method with Earth in circular orbit aroudn the Sun
       Outputs x,y,z,vx,vy,vz, for each planet. vx,vy,vz in AU/yr.
  -V   Velocity Verlet method with Earth in circular orbit aroudn the Sun
       Outputs x,y,z,vx,vy,vz, for each planet. vx,vy,vz in AU/yr.
  -T   Time Euler's and Vel. Verlet with preset n=10^1 to 10^9. Outputs
       result averaged over 10 runs in csv format: n,euler,verlet
  -C   Creator mode. Experiment with Sun-Earth-Jupiter system. Select Earth's
       velocity [AU/yr.] with -v, include and specify Jupiter mass with -j,
       play with the 1/r^b relationship between radius and force
       with -b.
       Outputs total kinetic,potential energy, angular momentum,
       x,y,z,vx,vy,vz,ep,ek,am, for each planet. vx,vy,vz in AU/yr.
  -P   Calculates the perihelion precession of Mercury. Include gen. rel. with -r.
  -n   Number of integration points per year.
  -y   Years to simulate
  -b   1/r^b relationship between radius and gravitational force (only with -C)
  -v   Earth's velocity vy (only with -C)
  -j   Enables Jupiter under -C and specifies its mass multiplier
  -r   Includes general relativistic correction to Newton's grav. law. (only with -P)

Example:
$ ./project3 -C -n 100 -y 2 -v 5
ek,ep,am,0x,0y,0z,0vx,0vy,0vz,0ep,0ek,0am,1x,1y,1z,1vx,1vy,1vz,1ep,1ek,1am,
3.7500000000000003e-05,-1.1843525281307230e-04,1.5000000000000000e-05,0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,-1.1843525281307230e-04,0.0000000000000000e+00,0.0000000000000000e+00,1.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,5.0000000000000000e+00,0.0000000000000000e+00,-1.1843525281307230e-04,3.7500000000000003e-05,1.5000000000000000e-05
...
```