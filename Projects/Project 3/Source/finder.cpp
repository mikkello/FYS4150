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
	cout << "Press 3 to use Velocity Verlet for Earth's escape velocity calculation." << endl;
	cout << "Press 4 to use Velocity Verlet for the Earth-Jupiter-Sun system." << endl;
	cout << "Press 5 to use Velocity Verlet for the full Solar system." << endl;
	cout << "Press 6 to use Velocity Verlet for the perihelion precession of Mercury." << endl;
	cout << "Press 7 to abort." << endl;
	cin >> integer;
	cout << "Positions will be printed to positions" << integer << ".xyz" << endl;
	cout << "Variables will be printed to variables.txt" << endl;

		
	clock_t start, finish;
	
	if (integer == 1) {
		cout << "You have chosen Euler Forward for the Earth-Sun system." << endl;
		cout << "Calculating orbits..." << endl;
		

		solarSystem.createCelestialBody(vec3(0, 0, 0), vec3(0, 0, 0), 1.0);
		solarSystem.createCelestialBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 3e-6);
		start = clock();
		Euler solvers(dt);
		solarSystem.writeToFile("positions1.xyz");

		for (int timestep = 0; timestep<numTimesteps; timestep++) {
			solvers.integrateOneStep(solarSystem);
			solarSystem.writeToFile("positions1.xyz");
		}
		
		
	}

	else if (integer == 2) {
		cout << "You have chosen Velocity Verlet for the Earth-Sun system." << endl;
		cout << "Calculating orbits..." << endl;

		solarSystem.createCelestialBody(vec3(0, 0, 0), vec3(0, 0, 0), 1.0);
		solarSystem.createCelestialBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 3e-6);
		start = clock();
		Verlet solvers(dt);
		solarSystem.writeToFile("positions2.xyz");

		for (int timestep = 0; timestep<numTimesteps; timestep++) {
			solvers.VelocityVerlet(solarSystem);
			solarSystem.writeToFile("positions2.xyz");
		}
		
	}

	else if (integer == 3) {
		cout << "You have chosen Velocity Verlet for Earth's escape velocity calculation." << endl;
		cout << "Calculating orbits..." << endl;

		double k = 1.47;

		solarSystem.createCelestialBody(vec3(0, 0, 0), vec3(0, 0, 0), 1.0);
		solarSystem.createCelestialBody(vec3(1, 0, 0), vec3(0, 2 * k * M_PI, 0), 3e-6);
		start = clock();
		Verlet solvers(dt);
		solarSystem.writeToFile("positions3_k_1.47.xyz");

		for (int timestep = 0; timestep<numTimesteps; timestep++) {
			solvers.VelocityVerlet(solarSystem);
			solarSystem.writeToFile("positions3_k_1.47.xyz");
		}

	}

	else if (integer == 4) {
		cout << "You have chosen Velocity Verlet for Earth-Jupiter-Sun system." << endl;
		cout << "Calculating orbits..." << endl;

		// Sun
		solarSystem.createCelestialBody(vec3(0, 0, 0), vec3(0, 0, 0), 1.0);
		// Earth
		solarSystem.createCelestialBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 3e-6);
		// Jupiter
		solarSystem.createCelestialBody(vec3(5.20, 0.0, 0.0), vec3(0.0, 0.87*M_PI, 0.0), 1000 * 0.95e-3);
		start = clock();
		Verlet solvers(dt);
		solarSystem.writeToFile("positions4_x1000.xyz");

	    for (int timestep = 0; timestep<numTimesteps; timestep++) {
			solvers.VelocityVerlet(solarSystem);
			solarSystem.writeToFile("positions4_x1000.xyz");
		}
	

	}

	else if (integer == 5) {
		cout << "You have chosen Velocity Verlet for the full Solar system." << endl;
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
		solarSystem.writeToFile("positions5.xyz");

		for (int timestep = 0; timestep<numTimesteps; timestep++) {
			solvers.VelocityVerlet(solarSystem);
			solarSystem.writeToFile("positions5.xyz");
		}


	}
	else if (integer == 6) {
		cout << "You have chosen Velocity Verlet for the perihelion precession of Mercury." << endl;
		cout << "Calculating orbits..." << endl;

		// Sun 
		solarSystem.createCelestialBody(-1.2E-7*vec3(0.3075, 0, 0), -1.2E-7*vec3(0, 12.44, 0), 1.0);
		// Mercury
		solarSystem.createCelestialBody(vec3(0.3075, 0, 0), vec3(0, 12.44, 0), 1.2E-7);
		start = clock();

		Verlet solvers(dt);
		ofstream Out;    
		Out.open("perihelion_GR_1e6.txt");

		//solarSystem.writeToFile("positions6.xyz");
		double x = solarSystem.bodies().at(0).position.x() - solarSystem.bodies().at(1).position.x();
		double y = solarSystem.bodies().at(0).position.y() - solarSystem.bodies().at(1).position.y();
		double r = solarSystem.bodies().at(0).position.length() - solarSystem.bodies().at(1).position.length();
		//cout << r << endl;

		double thetaPrev = 0; double theta = 0;
		double rPrevPrev = 0; double rPrev = 0.1;
		int i = 0;

		for (int timestep = 0; timestep<numTimesteps; timestep++) {
			solarSystem.writeToFile("positions6.xyz");
			solvers.Mercury_VelocityVerlet(solarSystem);
			x = solarSystem.bodies().at(0).position.x() - solarSystem.bodies().at(1).position.x();
			y = solarSystem.bodies().at(0).position.y() - solarSystem.bodies().at(1).position.y();
			theta = atan2(y, x);
			r = solarSystem.bodies().at(0).position.length() - solarSystem.bodies().at(1).position.length();

			if ((r > rPrev) && (rPrev < rPrevPrev)) {
				cout << "Perihelion angle = " << thetaPrev * 180 / M_PI * 3600 << endl;
				Out << setiosflags(ios::showpoint | ios::uppercase);
				Out << setprecision(10) << setw(20) << thetaPrev << endl;
				i = i + 1;
				thetaPrev = theta; rPrevPrev = rPrev; rPrev = r;

			}
			else {
				thetaPrev = theta; rPrevPrev = rPrev; rPrev = r;
			}

		}
		Out.close();
		cout << i << endl;


	}

	else if (integer == 7) {
		cout << "abort" << endl;
		terminate();
	}
	else {
		cout << "This number is not an option." << endl;
	}

	finish = clock();
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
	return Years;
}