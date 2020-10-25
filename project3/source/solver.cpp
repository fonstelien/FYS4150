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


/* Solve with given number of steps and time step dt. Trajectories and velocities are */
/* logged in mat flight_log. */
void Solver::solve(long long int steps, double dt) {
  vec r, v;
  double r3;
  Planet *p1, *p2;
  int num_planets = planets.size();
  
  flight_log = mat(steps+1, num_planets*9+3);  // all x,y,z,vx,vy,vz,ep,ek,am + tot. ek,ep,am

  // Log initial position
  for (int i = 0; i < num_planets; i++) {
    p1 = planets[i];
    flight_log(0, 3+9*i+0) = p1->pos(0);
    flight_log(0, 3+9*i+1) = p1->pos(1);
    flight_log(0, 3+9*i+2) = p1->pos(2);
    flight_log(0, 3+9*i+3) = p1->vel(0);
    flight_log(0, 3+9*i+4) = p1->vel(1);
    flight_log(0, 3+9*i+5) = p1->vel(2);
    flight_log(0, 3+9*i+6) = potential_energy(p1);
    flight_log(0, 3+9*i+7) = p1->kinetic_energy();
    flight_log(0, 3+9*i+8) = p1->angular_momentum();
  }
  flight_log(0, 0) = total_kinetic_energy();
  flight_log(0, 1) = total_potential_energy();
  flight_log(0, 2) = total_angular_momentum();
  
  // Initiate accelerations
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

  // Solve 
  for (long long int s = 1; s < steps+1; s++) {
    // Update position
    for (int i = 0; i < num_planets; i++) 
      planets[i]->update_pos(dt);

    // Abort if planet crashed
    for (int i = 0; i < num_planets; i++) {
      p1 = planets[i];
      if (p1->crashed) {
	cerr << "error: planet " << i << " crashed at " << s << " iterations." << endl;
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
      flight_log(s, 3+9*i+0) = p1->pos(0);
      flight_log(s, 3+9*i+1) = p1->pos(1);
      flight_log(s, 3+9*i+2) = p1->pos(2);
      flight_log(s, 3+9*i+3) = p1->vel(0);
      flight_log(s, 3+9*i+4) = p1->vel(1);
      flight_log(s, 3+9*i+5) = p1->vel(2);
      flight_log(s, 3+9*i+6) = potential_energy(p1);
      flight_log(s, 3+9*i+7) = p1->kinetic_energy();
      flight_log(s, 3+9*i+8) = p1->angular_momentum();
    }
    flight_log(s, 0) = total_kinetic_energy();
    flight_log(s, 1) = total_potential_energy();
    flight_log(s, 2) = total_angular_momentum();
    
  }
}


/* Returns csv format header string for the flight_log */
string Solver::csv_header() {
  string header = "ek,ep,am,";
  string s;
  for (int i = 0; i < (int) planets.size(); i++) {
    s = to_string(i);
    header += s + "x," + s + "y," + s + "z," + s + "vx," + s + "vy," + s + "vz,";
    header += s + "ep," + s + "ek," + s + "am,";
  }

  return header;
}


/* Total potential energy of Planet planet in sun-masses*AU^2/yr.^2 */
double Solver::potential_energy(Planet *planet) {
  Planet *other;
  double r;
  double ep = 0.;

  for (int i = 0; i < (int) planets.size(); i++) {
    other = planets[i];
    if (other == planet)
      continue;
    
    r = norm(planet->pos - other->pos);
    ep += -FPS*planet->m*other->m/(beta-1)/pow(r, beta-1);
  }

  return ep;
}


/* Total potential energy of all Planets in planets */
double Solver::total_potential_energy() {
  double ep = 0.;

  for (int i = 0; i < (int) planets.size(); i++)
    ep += potential_energy(planets[i]);

  return .5*ep;  // multiply with one half to avoid counting energy twice
}


/* Total kinetic energy of all Planets in planets */
double Solver::total_kinetic_energy() {
  double ek = 0.;

  for (int i = 0; i < (int) planets.size(); i++)
    ek += planets[i]->kinetic_energy();

  return ek;
}


/* Total angular momentum of all Planets in planets */
double Solver::total_angular_momentum() {
  double am = 0.;

  for (int i = 0; i < (int) planets.size(); i++)
    am += planets[i]->angular_momentum();

  return am;
}
