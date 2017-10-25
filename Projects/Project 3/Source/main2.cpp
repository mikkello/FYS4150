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
	
	int Steps_Year = Finder_Steps_Year();
	int Years = Finder_Years();
	double dt = 1.0 / float(Steps_Year);
	int numTimesteps = Years*Steps_Year;
	int num_planets = find_solver(dt, numTimesteps);
	

	//Saving variable to file
	ofstream variables;
	variables.open("variables1.txt");
	variables << setiosflags(ios::showpoint | ios::uppercase);
	variables << "N_timesteps" << " " << "N_planets" << " " << "N_years" << " " << "dt" << endl;
	variables << setprecision(10) << setw(20) << numTimesteps << " " << num_planets <<  " " << Years << " " <<  dt << endl;
	variables.close();

	

	return 0;
}


