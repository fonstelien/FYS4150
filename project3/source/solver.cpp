#include "utils.hpp"

using namespace arma;
using namespace std;


/* Adds Planet planet to the back of the planets collection */
void Solver::add(Planet *planet) {
  planets.push_back(planet);
}


/* Builds a solar system from mat system */
void Solver::build(mat system) {
  double mass;
  vec pos(3), vel(3);
  
  for (unsigned int i = 0; i < system.n_rows; i++) {
    mass = system(i,0);
    pos(0) = system(i,1);
    pos(1) = system(i,2);
    pos(2) = system(i,3);
    vel(0) = system(i,4);
    vel(1) = system(i,5);
    vel(2) = system(i,6);
    vel *= 365;
    add(new Planet(mass, pos, vel));
  }
}


/* Frees dynamically allocated memory pointed to by planets */
void Solver::nuke() {
  int num_planets = planets.size();
  for (int i = 0; i < num_planets; i++)
    delete planets[i];
  planets.clear();
}


/* Calculates the denominator in the gravitational force equation. */
/* This is where the varying beta and general relativistic model is implemented. */
double Solver::force_denominator(Planet *planet, vec r) {
  double rn = norm(r);
  double denom = pow(rn, beta+1);

  if (general_relativistic) {
    double ln = norm(cross(r, planet->vel));
    denom /= 1 + 3*ln*ln/(rn*rn*LIGHTPSPEED*LIGHTPSPEED);
  }
  
  return denom;
}


/* Solve with given number of steps and time step dt. Trajectories and velocities are */
/* logged in mat flight_log. */
void Solver::solve(int steps, double dt) {
  vec r, v;
  double r3;
  Planet *p1, *p2;
  int num_planets = planets.size();
  
  flight_log = mat(steps+1, num_planets*6);

  // Log initial position
  for (int i = 0; i < num_planets; i++) {
    p1 = planets[i];
    flight_log(0, 6*i+0) = p1->pos(0);
    flight_log(0, 6*i+1) = p1->pos(1);
    flight_log(0, 6*i+2) = p1->pos(2);
    flight_log(0, 6*i+3) = p1->vel(0);
    flight_log(0, 6*i+4) = p1->vel(1);
    flight_log(0, 6*i+5) = p1->vel(2);
  }

  // Initiate accelerations
  for (int i = 0; i < num_planets-1; i++) {
    p1 = planets[i];
    for (int j = i+1; j < num_planets; j++) {
      p2 = planets[j];
      r = p2->pos - p1->pos;
      r3 = force_denominator(p1, r);
      p1->update_acc(p2, r, r3);
      r3 = force_denominator(p2, r);
      p2->update_acc(p1, -r, r3);
    }
  }

  // Solve 
  for (int s = 1; s < steps+1; s++) {
    // Update position
    for (int i = 0; i < num_planets; i++) 
      planets[i]->update_pos(dt);

    // Abort if planet crashed
    for (int i = 0; i < num_planets; i++) {
      p1 = planets[i];
      if (p1->crashed) {
	cerr << "error: planet " << i << " crashed." << endl;
	flight_log = flight_log.rows(0,s);
	return;
      }
    }

    // Update accelerations at new position
    for (int i = 0; i < num_planets-1; i++) {
      p1 = planets[i];
      for (int j = i+1; j < num_planets; j++) {
	p2 = planets[j];
	r = p2->pos - p1->pos;
	r3 = pow(norm(r), beta+1);
	p1->update_acc(p2, r, r3);
	p2->update_acc(p1, -r, r3);
      }
    }

    // Update velocity
    for (int i = 0; i < num_planets; i++)
      planets[i]->update_vel(dt);

    // Log new position and velocity
    for (int i = 0; i < num_planets; i++) {
      p1 = planets[i];
      flight_log(s, 6*i+0) = p1->pos(0);
      flight_log(s, 6*i+1) = p1->pos(1);
      flight_log(s, 6*i+2) = p1->pos(2);
      flight_log(s, 6*i+3) = p1->vel(0);
      flight_log(s, 6*i+4) = p1->vel(1);
      flight_log(s, 6*i+5) = p1->vel(2);
    }
  }
}


/* Returns csv format header string for the flight_log */
string Solver::csv_header() {
  string header = "";
  string s;
  for (int i = 0; i < (int) planets.size(); i++) {
    s = to_string(i);
    header += s + "x," + s + "y," + s + "z," + s + "vx," + s + "vy," + s + "vz,";
  }

  return header;
}


/* Potential energy between Planets p1, p2 in sun-masses*AU^2/yr.^2 */
double Solver::potential_energy(Planet *p1, Planet *p2) {
  double r = norm(p1->pos - p2->pos);
  return -FPS*p1->m*p2->m/(beta-1)/pow(r, beta-1);
}
