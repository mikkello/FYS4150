#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "solarsystem.h"
#include "euler.h"
#include "verlet.h"
#include "finder.h"
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#define _USE_MATH_DEFINES

using namespace std;
int integer;

////////////////////////////////////////////////////////
// Function that asks user for which problem to solve //
int find_solver(double dt, int numTimesteps) {
	SolarSystem solarSystem;
	double AUday_to_AUyear = 365.2422; double kgtosolarmass = 5.02734e-31; // Unit conversion
	
	cout << "Which problem do you want to solve?" << endl;
	cout << "Press 1 to use Forward Euler for the Earth-Sun system." << endl;
	cout << "Press 2 to use Velocity Verlet for the Earth-Sun system." << endl;
	cout << "Press 3 to use Velocity Verlet for Earth's escape velocity calculation." << endl;
	cout << "Press 4 to use Velocity Verlet for the Earth-Jupiter-Sun system." << endl;
	cout << "Press 5 to use Velocity Verlet for the full Solar System." << endl;
	cout << "Press 6 to use Velocity Verlet for the perihelion precession of Mercury (Newtonian)." << endl;
	cout << "Press 7 to use Velocity Verlet for the perihelion precession of Mercury (GR)." << endl;
	cout << "Press 8 to abort." << endl;
	cin >> integer;
	cout << " " << endl;
	cout << "Positions will be written to:" << " " << "Positions_" << std::to_string(integer) << "_dt_" << std::to_string(dt) << ".xyz" << endl;
	cout << "Variables will be written to:" << " " << "Variables_" << std::to_string(integer) << "_dt_" << std::to_string(dt) << ".txt" <<endl;
	cout << " " << endl;
	
	
	clock_t start, finish;
	if (integer == 1) {
		cout << "You have chosen Euler Forward for the Earth-Sun system (" << integer << ")." << endl;
		cout << "Calculating orbits..." << endl;
		
		solarSystem.createCelestialBody(vec3(0, 0, 0), vec3(0, 0, 0), 1.0); // Sun
		solarSystem.createCelestialBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 3e-6); // Earth, circular orbit
		
		start = clock();
		Euler solvers(dt);
		solarSystem.writeToFile("Positions_" + std::to_string(integer) + "_dt_" + std::to_string(dt) + ".xyz");

		for (int timestep = 0; timestep<numTimesteps; timestep++) {
			solvers.integrateOneStep(solarSystem);
			solarSystem.writeToFile("Positions_" + std::to_string(integer) + "_dt_" + std::to_string(dt) + ".xyz");
		}
		
		
	}

	else if (integer == 2) {
		cout << "You have chosen Velocity Verlet for the Earth-Sun system (" << integer << ")." << endl;
		cout << "Calculating orbits..." << endl;

		solarSystem.createCelestialBody(vec3(0, 0, 0), vec3(0, 0, 0), 1.0); // Sun
		solarSystem.createCelestialBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 3e-6); // Earth, circular orbit
		
		start = clock();
		Verlet solvers(dt);
		solarSystem.writeToFile("Positions_" + std::to_string(integer) + "_dt_" + std::to_string(dt) + ".xyz");

		for (int timestep = 0; timestep<numTimesteps; timestep++) {
			solvers.VelocityVerlet(solarSystem);
			solarSystem.writeToFile("Positions_" + std::to_string(integer) + "_dt_" + std::to_string(dt) + ".xyz");
		}
	}

	else if (integer == 3) {
		cout << "You have chosen Velocity Verlet for Earth's escape velocity calculation (" << integer << ")." << endl;
		cout << "Calculating orbits..." << endl;

		double k = 1.41; // v0 = k * 2 * pi // Initial velocity constant

		solarSystem.createCelestialBody(vec3(0, 0, 0), vec3(0, 0, 0), 1.0); // Sun
		solarSystem.createCelestialBody(vec3(1, 0, 0), vec3(0, 2 * k * M_PI, 0), 3e-6); // Earth, escaping

		start = clock();
		Verlet solvers(dt);
		solarSystem.writeToFile("Positions_" + std::to_string(integer) + "_k_" + std::to_string(k) + "_dt_" + std::to_string(dt) + ".xyz");

		for (int timestep = 0; timestep<numTimesteps; timestep++) {
			solvers.VelocityVerlet(solarSystem);
			solarSystem.writeToFile("Positions_" + std::to_string(integer) + "_k_" + std::to_string(k) + "_dt_" + std::to_string(dt) + ".xyz");
		}
	}

	else if (integer == 4) {
		cout << "You have chosen Velocity Verlet for Earth-Jupiter-Sun system (" << integer << ")." << endl;
		cout << "Calculating orbits..." << endl;
		int jupiter_mass_const = 1; // Change Jupiter's mass by multiplying with this constant

		solarSystem.createCelestialBody(vec3(0, 0, 0), vec3(0, 0, 0), 1.0); // Sun
		solarSystem.createCelestialBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 3e-6); // Earth 
		solarSystem.createCelestialBody(vec3(5.20, 0.0, 0.0), vec3(0.0, 0.87*M_PI, 0.0), jupiter_mass_const * 0.95e-3); // Jupiter
		
		start = clock();
		Verlet solvers(dt);
		solarSystem.writeToFile("Positions_" + std::to_string(integer) + "_mX" + std::to_string(jupiter_mass_const) + "_dt_" + std::to_string(dt) + ".xyz");

	    for (int timestep = 0; timestep<numTimesteps; timestep++) {
			solvers.VelocityVerlet(solarSystem);
			solarSystem.writeToFile("Positions_" + std::to_string(integer) + "_mX" + std::to_string(jupiter_mass_const) + "_dt_" + std::to_string(dt) + ".xyz");
		}
	}

	else if (integer == 5) {
		cout << "You have chosen Velocity Verlet for the full Solar System (" << integer << ")." << endl;
		cout << "Calculating orbits..." << endl;

		// Sun
		solarSystem.createCelestialBody(vec3(0, 0, 0), vec3(-2.013685881178907E-06, 6.833576056500457E-06, 4.125365253557130E-08)*AUday_to_AUyear, 1.0);
		// Mercury
		solarSystem.createCelestialBody(vec3(-3.767421263697904E-01,  2.668580947011324E-02,  3.663085930504686E-02),
		                                     vec3(-7.569397216231175E-03, -2.686580344697246E-02, -1.501528485574016E-03)*AUday_to_AUyear,
		                                    1.2e-7);
		// Venus
		solarSystem.createCelestialBody(vec3(2.392368478680558E-01, -6.847999876013139E-01, -2.319637912277905E-02),
		                                     vec3(1.899757861343120E-02,  6.491099711740283E-03, -1.007444652406855E-03)*AUday_to_AUyear,
		                                     2.45e-6);
		// Earth
		solarSystem.createCelestialBody(vec3(9.152108693600010E-01,  4.058839764119083E-01, -1.771902422530673E-04),
		                                    vec3(-7.226100108439347E-03,  1.567554554659565E-02, -1.035367359990257E-07)*AUday_to_AUyear,
		                                     3e-6);
		// Mars
		solarSystem.createCelestialBody(vec3(1.178122338122247E+00, -7.242078750401497E-01, -4.423390202799733E-02),
		                                     vec3(7.899668162477799E-03,  1.310156103432712E-02,  8.052516290965436E-05)*AUday_to_AUyear,
		                                     3.3e-7);
		// Jupiter
		solarSystem.createCelestialBody(vec3(-5.427455884765774E+00, -4.678991629217762E-01,  1.233233582953724E-01),
		                                    vec3(5.597320182818523E-04, -7.162261501832823E-03,  1.717650487472710E-05)*AUday_to_AUyear,
		                                    9.95e-4);
		// Saturn
		solarSystem.createCelestialBody(vec3(-2.256792396133917E+00, -9.777367032001933E+00,  2.598201193469362E-01),
		                                     vec3(5.130301737533382E-03, -1.272889194718095E-03, -1.817149598022347E-04)*AUday_to_AUyear,
		                                     2.75e-4);
		// Uranus
		solarSystem.createCelestialBody(vec3(1.846018493184810E+01,  7.568337250358790E+00, -2.110464444997344E-01),
		                                     vec3(-1.520797540482449E-03,  3.455811674219709E-03,  3.259634893752387E-05)*AUday_to_AUyear,
		                                     4.4e-5);
		// Neptune
		solarSystem.createCelestialBody(vec3(2.826297032275447E+01, -9.916332407218558E+00, -4.471405099000715E-01),
		                                     vec3(1.018822778770138E-03,  2.980778775480523E-03, -8.516397608222377E-05)*AUday_to_AUyear,
		                                     5.15e-5);

		start = clock();
		Verlet solvers(dt);
		solarSystem.writeToFile("Positions_" + std::to_string(integer) + "_dt_" + std::to_string(dt) + ".xyz");

		for (int timestep = 0; timestep<numTimesteps; timestep++) {
			solvers.VelocityVerlet(solarSystem);
			solarSystem.writeToFile("Positions_" + std::to_string(integer) + "_dt_" + std::to_string(dt) + ".xyz");
		}
	}

	else if (integer == 6) {
		cout << "You have chosen Velocity Verlet for the perihelion precession of Mercury (" << integer << ")." << endl;
		cout << "Calculating orbits..." << endl;

		solarSystem.createCelestialBody(vec3(0, 0, 0), vec3(0, 0, 0), 1); // Sun
		solarSystem.createCelestialBody(vec3(0.3075, 0, 0), vec3(0, 12.44, 0), 1.2e-7); // Mercury

		start = clock();
		double time = 0;
		Verlet solvers(dt);
		
		ofstream Outfile;
		Outfile.open("Perihelion_" + std::to_string(integer) + "_dt_" + std::to_string(dt) + ".txt");
		Outfile << "Perihelion angle (arcseconds)" << " " << "Time" << endl;

		for (int i = 0; i < (numTimesteps) / 2; i++) {
			double pos_min2 = solarSystem.bodies()[1].position.length();
			solvers.VelocityVerlet(solarSystem);
			time += dt;

			double pos_min1 = solarSystem.bodies()[1].position.length();
			vec3 periPos = solarSystem.bodies()[1].position;

			solvers.VelocityVerlet(solarSystem);
			time += dt;

			double pos = solarSystem.bodies()[1].position.length();

			// Writing arcseconds and time to file if located at the perihelion
			if ((pos_min1 < pos_min2) && (pos_min1 < pos)) {
				double tanTheta = periPos[1] / periPos[0];
				double theta = atan(tanTheta);
				double arcsec = (theta * 3600 * 180) / M_PI; // Converting radians to arcseconds

				Outfile << arcsec << " " << (time - dt) << endl;
			}
		}
		Outfile.close();


	}
	else if (integer == 7) {
		cout << "You have chosen Velocity Verlet for the perihelion precession of Mercury (GR) (" << integer << ")." << endl;
		cout << "Calculating orbits..." << endl;
		 
		solarSystem.createCelestialBody(vec3(0, 0, 0), vec3(0, 0, 0), 1); // Sun
		solarSystem.createCelestialBody(vec3(0.3075, 0, 0), vec3(0, 12.44, 0), 1.2E-7); // Mercury

		start = clock();

		Verlet solvers(dt);
		ofstream Outfile;

		double time = 0;
		Outfile.open("Perihelion_" + std::to_string(integer) + "_dt_" + std::to_string(dt) + ".txt");
		Outfile << "Perihelion angle (arcseconds)" << " " << "Time" << endl;
		

		Verlet integrator(dt);
		for (int i = 0; i < (numTimesteps / 2) ; i++) {

			double pos_min2 = solarSystem.bodies()[1].position.length();
			integrator.relativistic_VelocityVerlet(solarSystem);
			time += dt;

			double pos_min1 = solarSystem.bodies()[1].position.length();
			vec3 periPos = solarSystem.bodies()[1].position;

			integrator.relativistic_VelocityVerlet(solarSystem);
			time += dt;
			double pos = solarSystem.bodies()[1].position.length();

			// Writing arcseconds and time to file if located at the perihelion
			if ((pos_min1 < pos_min2) && (pos_min1 < pos)) {

				double tanTheta = periPos[1] / periPos[0];
				double theta = atan(tanTheta);
				double arcsec = (theta * 3600 * 180) / M_PI; // Converting radians to arcseconds

				Outfile << arcsec << " " << (time - dt) << endl;
			}
		}
		Outfile.close();
	}
	
	else if (integer == 8) {
		cout << "abort" << endl;
		std::abort();
	}
	else {
		cout << "This number is not an option." << endl;
	}

	finish = clock(); // Timing
	double time = (finish - start); 

	cout << "Orbit calculation done." << endl;
	cout << "CPU time:" << time / CLOCKS_PER_SEC << " s." << endl;

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
	cout << " " << endl;
	return Years;
}

int Problem_Chosen() {
	// Returns the number according to which problem that is chosen to solve. 
	return integer;
}