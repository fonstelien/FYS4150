#include "utils.hpp"

using namespace arma;
using namespace std;

// Speed of light in AU per year
#define LIGHTPSPEED (2.99792458E+8 / 1.495978707E+11 * 60*60*24*365)

// One arc second
#define ARCSECOND (M_PI/180/60/60)

/* Runs Velocity Verlet with earth in circular orbit and logs */
/* trajectory and velocity in flight_log. */
double perihelion_of_mercury(double h, long long int n, bool general_relativistic) {
  vec pos(2), vel(2), acc0(2), acc(2);
  double rn, ln, gr;
  double x0, x1, xp, y0, y1, yp;
  double r0, r1, r2;
  double arcsecs;
  double c2 = LIGHTPSPEED*LIGHTPSPEED;
  
  pos = {0.3075,0.,0.};
  vel = {0.,12.44,0.};

  x0 = x1 = xp = pos(0);
  y0 = y1 = yp = pos(1);
  r0 = x0*x0 + y0*y0;
  r1 = r2 = -1.;  

  // Initializing
  rn = norm(pos);
  gr = 1.;
  if (general_relativistic) {
    ln = vel(0)*pos(1) - vel(1)*pos(0);
    gr = 1. + 3*ln*ln/(rn*rn*c2);
  }
  acc = -FPS/(rn*rn*rn)*gr*pos;
  
  for (long long int i = 0; i < n; i++) {
    // Update position and velocity
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

    // Find the perihelion
    x0 = pos(0);
    y0 = pos(1);
    r2 = r1;
    r1 = r0;
    r0 = x0*x0 + y0*y0;
    if (r1 < r0 && r1 < r2) {  // Was the former point the perihelion?
      xp = x1;
      yp = y1;
    }
    x1 = x0;
    y1 = y0;
  }
  
  arcsecs = atan(yp/xp)/ARCSECOND;
  return arcsecs;
}
