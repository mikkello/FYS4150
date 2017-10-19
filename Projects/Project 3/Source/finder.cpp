#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "solarsystem.h"
#include "euler.h"
#include "verlet.h"
#include <fstream>
#include <iomanip>
#include <cstdlib>
#define _USE_MATH_DEFINES

using namespace std;

// Function that asks user for which problem to solve
int find_solver(double dt, int numTimesteps) {
	SolarSystem solarSystem;
	double AUday_to_AUyear = 365.24212;
	double kgtosolarmass = 5.02734e-31;

	int integer;
	cout << "Which problem do you want to solve?" << endl;
	cout << "Press 1 to use Forward Euler for the Earth-Sun system." << endl;
	cout << "Press 2 to use Velocity Verlet for the Earth-Sun system." << endl;
	cout << "Press 3 to abort." << endl;
	cin >> integer;
	cout << "Positions will be printed to positions" << integer << ".xyz" << endl;
	cout << "Variables will be printed to variables" << integer << ".txt" << endl;

		
		
	if (integer == 1) {
		cout << "You have choosen Euler Forward for the Earth-Sun system." << endl;
		cout << "Calculating orbits..." << endl;
		solarSystem.createCelestialBody(vec3(0, 0, 0), vec3(0, 0, 0), 1.0);
		solarSystem.createCelestialBody(vec3(9.890331046925951E-01, 1.768079890757788E-01, -1.738715302893284E-04), 
vec3(-3.268395786841218E-03, 1.689265025904021E-02, -9.889230545039174E-07)*AUday_to_AUyear, 5.97219e24*kgtosolarmass);

		Euler solvers(dt);
		solarSystem.writeToFile("positions1.xyz");

		for (int timestep = 0; timestep<numTimesteps; timestep++) {
			solvers.integrateOneStep(solarSystem);
			solarSystem.writeToFile("positions1.xyz");
		}
		
	}

	else if (integer == 2) {
		cout << "You have choosen Euler Forward for the Earth-Sun system." << endl;
		cout << "Calculating orbits..." << endl;

		solarSystem.createCelestialBody(vec3(0, 0, 0), vec3(0, 0, 0), 1.0);
		solarSystem.createCelestialBody(vec3(9.890331046925951E-01, 1.768079890757788E-01, -1.738715302893284E-04), 
vec3(-3.268395786841218E-03, 1.689265025904021E-02, -9.889230545039174E-07)*AUday_to_AUyear, 5.97219e24*kgtosolarmass);

		Verlet solvers(dt);
		solarSystem.writeToFile("positions2.xyz");

		for (int timestep = 0; timestep<numTimesteps; timestep++) {
			solvers.VelocityVerlet(solarSystem);
			solarSystem.writeToFile("positions2.xyz");
		}

	}

	else if (integer == 3) {
		cout << "abort" << endl;
		terminate();
	}
	else {
		cout << "This number is not an option." << endl;
	}

	
	return solarSystem.bodies().size();
}

int Finder_Steps_Year() {
	// Finds time steps per year
	int Steps_Year;
	cout << "Time steps per year?" << endl;
	cin >> Steps_Year;
	return Steps_Year;
}
int Finder_Years() {
	// Finds years for the calculations
	int Years;
	cout << "How many years would you like to simulate?" << endl;
	cin >> Years;
	return Years;
}