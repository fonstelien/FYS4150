#include <time.h>
#include "utils.hpp"

using namespace arma;
using namespace std;

// One arc second
#define ARCSECOND (M_PI/180/60/60)

/* Runs Velocity Verlet with earth in circular orbit and logs */
/* trajectory and velocity in flight_log. */
double perihelion_of_mercury(double h, long long int n, bool general_relativistic) {
  long long int batch = 100000;
  int idx = -1;
  double rn, ln, gr;
  double c2 = LIGHTPSPEED*LIGHTPSPEED;
  mat trajectory(batch+1,2);
  mat section;
  vec pos(2), vel(2), acc0(2), acc(2);
  vec x, y, r;
  double arcsecs;

  if (batch > n)
    batch = n;
  
  pos = {0.3075,0.,0.};
  vel = {0.,12.44,0.};

  for (long long int b = 0; b < n; b += batch) {
    if (b+batch > n)
      batch = n - b;
    
    trajectory(0,0) = pos(0);
    trajectory(0,1) = pos(1);

    // Initializing
    rn = norm(pos);
    gr = 1.;
    if (general_relativistic) {
      ln = vel(0)*pos(1) - vel(1)*pos(0);
      gr = 1. + 3*ln*ln/(rn*rn*c2);
    }
    acc = -FPS/(rn*rn*rn)*gr*pos;

    for (int i = 1; i < batch+1; i++) {
      pos += h*vel + h*h/2*acc;
      acc0 = acc;

      rn = norm(pos);
      gr = 1.;
      if (general_relativistic) {
	ln = vel(0)*pos(1) - vel(1)*pos(0);
	gr = 1. + 3*ln*ln/(rn*rn*c2);
      }
      acc = -FPS/(rn*rn*rn)*gr*pos;

      vel += h/2*(acc + acc0);
    
      trajectory(i,0) = pos(0);
      trajectory(i,1) = pos(1);
    }

    if(trajectory(0,1) < 0 && trajectory(batch,1) > 0)
      section = trajectory;
  }
  
  x = section.col(0);
  y = section.col(1);
  r = square(x) + square(y);
  for (idx = r.n_elem-2; idx > 1; idx--)
    if (r[idx] < r[idx+1] && r[idx] < r[idx-1])
      break;
  
  arcsecs = atan(y[idx]/x[idx])/ARCSECOND;

  return arcsecs;
}
