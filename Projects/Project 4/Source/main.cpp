#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <string>

using namespace std;

double calculateE(int size, double J, double** spinMatrix);
double calcdE(double **spinMatrix, int i, int j, int L);
double calculateM(int L, double **spinMatrix);
void metropolisAlg(int numprocs, double T, double& E, double& M, double& E2, double& M2, double& cv, double& chi, int L, int N, double J, double** spinMatrix, string filename);
void printData(double L, double T, double N, double energy, double E2, double M2, double magnetization, double cv, double chi);
void normalizing(int numprocs, double T, int L, int N, double& E, double& E2, double& M, double& M2, double& cv, double& chi);

double** generateOrderedMat(int L);
double** generateRandMat(int L);



int main(int nargs, char* args[]) {
	int N = 1e5; 		// Monte Carlo Cycles
	int L = 100;			// Lattice size, LxL 
	double J = -1.0;	// Coupling constant
	double mintemp = 2.25; // Start temp
	double maxtemp = 2.35; // Stop temp
	double tempstep = 0.05; // Temp step
	 

	// Initializing matrix setup
	// Change to "ordered" for ordered lattice
	double** spins;
	spins = generateRandMat(L);
    string initialSpinState = "random";
	//spins = generateOrderedMat(L);
	//string initialSpinState = "ordered";

	// Paralellization:
	int numprocs, my_rank;
	double start, end;
	MPI_Init(&nargs, &args);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// Seed for RNG
	srand(time(NULL) + my_rank);

	// File for final values
	ofstream myfile;
	myfile.open("L" + std::to_string(L)+ "_" + std::to_string(N) + ".dat");

	// Main algorithm
	myfile << "E_tot" << " " << "total_M" << " " << "cv" << " " << "chi" << " " << "temp" << endl;
	for (double temp = mintemp; temp <= maxtemp; temp += tempstep) {

		// Declaring variables
		double E, M, E2, M2, cv, chi;
		double E_tot = 0; double total_M = 0;
		double E2_tot = 0; double M2_tot = 0;
		//double acc = 0;
		char filename[10000];
		sprintf(filename, "E_M_%s_T%.2f.dat", initialSpinState.c_str(), temp);

		MPI_Barrier(MPI_COMM_WORLD);
		start = MPI_Wtime();

		metropolisAlg(numprocs, temp, E, M, E2, M2, cv, chi, L, N, J, spins, string(filename)); // Initializing algorithm

		MPI_Reduce(&E, &E_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&M, &total_M, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&E2, &E2_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&M2, &M2_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		MPI_Barrier(MPI_COMM_WORLD); 
		end = MPI_Wtime();

		if (my_rank == 0) {
			normalizing(numprocs, temp, L, N, E_tot, E2_tot, total_M, M2_tot, cv, chi);
			printData(L, temp, N*numprocs, E_tot, E2_tot, M2_tot, total_M, cv, chi);
			myfile << E_tot << " " << total_M << " " << cv << " " << chi << " " << temp << " " << endl;
			cout << "CPU time: " << (end - start) << endl;
		}
		
	}
	myfile.close(); 
	MPI_Finalize();
	
	return 0;
}

double calculateM(int L, double **spinMatrix) {
	double magnetization = 0;
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			magnetization += spinMatrix[i][j];
		}
	}
	return fabs(magnetization);
}

double calculateE(int size, double J, double **spinMatrix) {
	double mean_E = 0;
	int L = size;
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {

			int up = j + 1; int down = j - 1;
			int left = i - 1; int right = i + 1;

			if (i == L - 1) right = 0;
			if (i == 0)   left = L - 1;
			if (j == L - 1) up = 0;
			if (j == 0)   down = L - 1;

			mean_E += (spinMatrix[i][j] * spinMatrix[i][up]) + (spinMatrix[i][j] * spinMatrix[i][down])
				+ (spinMatrix[i][j] * spinMatrix[left][j]) + (spinMatrix[i][j] * spinMatrix[right][j]);
		}
	}
	return J*mean_E*0.5;
}


double** generateOrderedMat(int L) {
	double **spins = new double*[L];
	for (int i = 0; i < L; i++) {
		spins[i] = new double[L];
	}
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			spins[i][j] = 1;
		}
	}
	return spins;
}

double** generateRandMat(int L) {
	double **spins = new double*[L];
	for (int i = 0; i < L; i++) {
		spins[i] = new double[L];
	}
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			int x = rand() % 2;
			if (x == 1) {
				spins[i][j] = 1;
			}
			else {
				spins[i][j] = -1;
			}
		}
	}
	return spins;
}


double calcdE(double **spinMatrix, int i, int j, int L) {
	int up = j + 1; int down = j - 1;
	int left = i - 1; int right = i + 1;
	if (j == L - 1) up = 0;
	if (j == 0)   down = L - 1;
	if (i == L - 1) right = 0;
	if (i == 0)   left = L - 1;

	double dE = 2 * (spinMatrix[i][j] * (spinMatrix[i][up] + spinMatrix[i][down] + spinMatrix[left][j] + spinMatrix[right][j]));
	return dE;
}

void metropolisAlg(int numprocs, double T, double& E, double& M, double& E2, double& M2, double& cv, double& chi, int L, int N, double J, double** spinMatrix, string filename) {
	// Uncomment to write to file for every cycle
	//ofstream myfile2;
	//myfile2.open(filename.c_str());

	// Calculating energy
	double energy = calculateE(L, J, spinMatrix);
	double magnetization = calculateM(L, spinMatrix); 
	double magneSum = 0; double energySum = 0;
	double energySqrSum = 0; double magneSqrSum = 0;
	int n = 0;
	double beta = 1./ T;
	int accept = 0;
	while (n < N) {
		for (int k = 0; k < L*L; k++) {
			// Picks a random atom
			int i = rand() % L; int j = rand() % L;
			double dE = calcdE(spinMatrix, i, j, L);
			magnetization = calculateM(L, spinMatrix);
			if (dE <= 0) {
				// Accepting new energy
				energy += dE;
				magnetization = calculateM(L, spinMatrix);
				accept += 1;
				// Flipping spin
				spinMatrix[i][j] *= (-1);
			}
			else {
				double w = exp(-beta*dE);
				double max = RAND_MAX;
				double r = rand() / ((double)max);
				if (r < w) {
					energy += dE;
					magnetization = calculateM(L, spinMatrix);
					accept += 1;
					// Flipping spin
					spinMatrix[i][j] *= (-1);
				}
			}
		}
		energySum += energy;
		energySqrSum += pow(energy, 2);
		magneSum += magnetization;
		magneSqrSum += pow(magnetization, 2);
		n += 1;
		
		//myfile2 << n << endl;
		//myfile2 << accept << endl;
		//myfile2 << energySum/((double)n*L*L) << " " << abs(magneSum/((double) n*L*L)) << endl;
		//myfile2 << energy << " " << magnetization << endl;

	}
	E = energySum; M = magneSum;
	E2 = energySqrSum; M2 = magneSqrSum;
	//acc = accept;

	//myfile2.close();		
}

void normalizing(int numprocs, double T, int L, int N, double& E, double& E2, double& M, double& M2, double& cv, double& chi) {
	double LSqr = pow(L, 2);
	E = E / (numprocs*N);
	E2 = E2 / (numprocs*N);
	M = M / (numprocs*N);
	M2 = M2 / (numprocs*N);
	cv = (E2 - pow(E, 2)) / pow(T, 2);
	cv /= LSqr;
	chi = (M2 - pow(M, 2)) / T;
	chi /= LSqr;
	E /= LSqr;
	E2 /= LSqr;
	M /= LSqr;
	M2 /= LSqr;
}

void printData(double L, double T, double N, double energy, double E2, double M2, double magnetization, double cv, double chi) {
	cout << endl;
	cout << "Loops: " << N << " Temperature: " << T << " L = " << L << endl;
	cout << setprecision(4) << "  Mean energy:               " << energy << endl;
	cout << setprecision(4) << "  Mean magnetization:        " << magnetization << endl;
	cout << setprecision(4) << "  Heat capacity:             " << cv << endl;
	cout << setprecision(4) << "  Susceptibility             " << chi << endl;
}

