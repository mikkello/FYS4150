// This project is a part of the course FYS4150 - Computational Physics @ UiO. 
// The code is based on an example by Anders Hafreager
// https://github.com/andeplane/solar-system

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include "solarsystem.h"
#include "euler.h"
#include "finder.h"
#include <cstdio>
using namespace std;

int main()
{
	// Timing
	clock_t start, finish;
	start = clock();

	int Steps_Year = Finder_Steps_Year();
	int Years = Finder_Years();
	double dt = 1.0 / float(Steps_Year);
	int numTimesteps = Years*Steps_Year;
	int num_planets = find_solver(dt, numTimesteps);

	// Saving variable to file
	ofstream variables;
	variables.open("variables.txt");
	variables << setiosflags(ios::showpoint | ios::uppercase);
	variables << setprecision(10) << setw(20) << numTimesteps << " " << num_planets << endl;
	variables.close();

	// Stop timing
	finish = clock();
	double time = (finish - start);
	
	cout << "Orbit calculation done." << endl;
	cout << "CPU time:" << time / CLOCKS_PER_SEC << " s." << endl;
	cout << "Positions written to file positions.xyz." << endl;
	cout << "Positions written to file variables.txt." << endl;
	

	return 0;
}


