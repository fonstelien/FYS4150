#include "utils.hpp"

using namespace arma;
using namespace std;

void Solver::add(Planet planet) {
  planets.push_back(planet);
}

void Solver::solve(int steps, double dt) {
  vec r;
  double r3;
  Planet *p1, *p2;
  int num_planets = planets.size();
  
  flight_log = mat(steps, num_planets*6);

  for (int s = 0; s < steps; s++) {
    for (int i = 0; i < num_planets-1; i++) {
      p1 = &planets[i];
      for (int j = i+1; j < num_planets; j++) {
	p2 = &planets[j];
	r = p2->pos - p1->pos;
	r3 = norm(r);
	r3 = r3*r3*r3;
	p1->update_acc0(p2, r, r3);
	p2->update_acc0(p1, -r, r3);
      }
    }

    for (int i = 0; i < num_planets; i++)
      planets[i].update_pos(dt);

    for (int i = 0; i < num_planets-1; i++) {
      p1 = &planets[i];
      for (int j = i+1; j < num_planets; j++) {
	p2 = &planets[j];
	r = p2->pos - p1->pos;
	r3 = norm(r);
	r3 = r3*r3*r3;
	p1->update_acc1(p2, r, r3);
	p2->update_acc1(p1, -r, r3);
      }
    }
  
    for (int i = 0; i < num_planets; i++)
      planets[i].update_vel(dt);

    for (int i = 0; i < num_planets; i++) {
      p1 = &planets[i];
      flight_log(s, 6*i+0) = p1->pos(0);
      flight_log(s, 6*i+1) = p1->pos(1);
      flight_log(s, 6*i+2) = p1->pos(2);
      flight_log(s, 6*i+3) = p1->vel(0);
      flight_log(s, 6*i+4) = p1->vel(1);
      flight_log(s, 6*i+5) = p1->vel(2);
    }
  }
}


void Solver::to_csv() {
  for (int i = 0; i < (int) planets.size(); i++)
    cout << i << "x," << i << "y," << i << "z," << i << "vx," << i << "vy," << i << "vz,";
  cout << endl;
  flight_log.save(cout, csv_ascii);
}
