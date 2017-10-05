//Eigenvalue solver - Jacobi's algorithm.
//In this project, eigenvalues and eigenvectors for electrons trapped in harmonic oscilltor potentials are found.
//The progject is part of FYS4150 - Computational Physics @ UiO.
//The coded algorithm is inspiried by the lecture notes of Morten Hjorth-Jensen
//Mikkel Christensen - October 2017.

#include <armadillo>
#include <iostream>
#include <ctime>
#include <iomanip>

using namespace std;
using namespace arma;

///////////////////////
//Potential functions//
double potential_one_e(double rho){
	return pow(rho, 2);
}

double potential_two_e_no_interaction(double rho, double omega){
	return pow(omega, 2)*pow(rho, 2);
}

double potential_two_e_interaction(double rho, double omega){
	return pow(omega, 2)*pow(rho, 2) + 1.0 / rho;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Jacobi rotation algorithm. Inputs: matrix A, eigenvector matrix R, row k, column l, number of iterations N.//
void Jac_Rot(mat &A, mat &R, int k, int l, int N)
{
	double s, c; // s = sine(theta), c = cos(theta)
	double t, tau; // t = tan(theta), tau = cos(2theta)
// Finding theta 
	if (A(k, l) != 0)
	{
		tau = (A(l, l) - A(k, k)) / (2 * A(k, l));
		if (tau >= 0)
		{
			t = 1.0 / (tau + sqrt(1.0 + pow(tau, 2)));
		}
		else
		{
			t = -1.0 / (-tau + sqrt(1 + pow(tau, 2)));
		}
		c = 1.0 / sqrt(1 + pow(t, 2));
		s = c*t;
	}
	else
	{
		c = 1.0;
		s = 0.0;
	}

// Rotation
	double Akk = A(k, k);
	A(k, k) = A(k, k)*pow(c, 2) - 2 * A(k, l)*c*s + A(l, l)*pow(s, 2);
	A(l, l) = A(l, l)*pow(c, 2) + 2 * A(k, l)*c*s + Akk*pow(s, 2);
	A(k, l) = A(l, k) = 0;
	for (int i = 0; i < N - 1; i++) 
	{
		if (i != k && i != l) {
			double Aik = A(i, k);
			A(i, k) = c*A(i, k) - s*A(i, l);
			A(i, l) = c*A(i, l) + s*Aik;
			A(k, i) = A(i, k);
			A(l, i) = A(i, l);
		}
// Eigenvector matrix
		double R_il = R(i, l);
		double R_ik = R(i, k);
		R(i, l) = c*R_il + s*R_ik;
		R(i, k) = -s*R_il + c*R_ik;
	}
}
// End of rotation algorithm 



//////////////////////////
//Largest element finder//
void largest_element(mat A, int* kmax, int* lmax) {
	int N;
// Quadratic matrix test
	if (A.n_rows!= A.n_cols) {           
		throw invalid_argument("Error: N_rows does not equal N_cols.");
	}
	else {
		N = A.n_rows;
	}
	double max = 0;
	for (int k = 0; k < N; k++) {
		for (int l = 0; l < N; l++) {
			if (k != l) {
				if (abs(A(k, l)) >= max) {
					*kmax = k; *lmax = l;
					max = abs(A(k, l));
				}
			}
		}
	}
}
// End of largest element finder



//////////////////////
//Orthogonality test//
void ortho_test(mat R, double eps) {
	int R_cols = R.n_cols;
	int R_rows = R.n_rows;
// Column vectors
	vec Col_1 = zeros<vec>(R_rows);
	vec Col_2 = zeros<vec>(R_rows);
	for (int c_1 = 0; c_1 < R_cols; c_1++) {
		for (int i = 0; i < R_rows; i++) {
			// Column i
			Col_1(i) = R(c_1, i); 
		}
		for (int c_2 = 0; c_2 < R_cols; c_2++) {
			if (c_1 != c_2) {
				for (int j = 0; j < R_rows; j++) {
					// Column j
					Col_2(j) = R(c_2, j); 
				}
				double dotproduct = dot(Col_1, Col_2); // Dot product column i and j
				if (dotproduct > eps) {
					cout << "Dot product is not conserved." << endl;
				}
			}
		}
	}
}
// End of orthogonality test



////////////////////////
//Laegest element test//
void test_largest_element(int k, int l) {
	int N;
//Test if element is on the diagonal
	if (k == l) {
		throw  invalid_argument("k = l. Largest element is on the diagonal.");
	}
	else if (k >= l) {
		N = k;
	}
	else {
		N = l;
	}
	mat A = zeros<mat>(N + 1, N + 1);
// Place 1 at (k,l)
	A(k, l) = 1;
// Store previous k and l
	int k_o = k; int l_o = l;
// Test the new matrix
	largest_element(A, &k, &l);
// Test if the functions finds the correct location 
	if (k != k_o && l != l_o) {
		throw invalid_argument("Largest_element function failed");
	}
}
// End of largest element test




////////////////////////////////////////////////////////////////////////////////////////////////
//Finding eigenvalues and eigenvectors for electrons trapped in Harmonic Oscilaltor potensials//
int main() {

	int N = 10; // Grid points. Varied in project.
	clock_t start_1, finish_1, start_2, finish_2;
	double* rho = new double[N + 1];
	rho[0] = 0.0;
	rho[N] = 10.0; // Rho_max. Varied in project.
	double h = (rho[N] - rho[0]) / N; // Step length
	double omega = 0.01; // Varied in project: [0.01, 0.5, 1, 5]

// Rho
	for (int i = 1; i < N; i++)
	{
		rho[i] = rho[0] + i*h;
	}


// Potential matrix A
	mat A = zeros<mat>(N - 1, N - 1);

	for (int i = 0; i < N - 1; i++)
	{
		// Change potential here!
		A(i, i) = (2.0 / pow(h, 2)) + potential_two_e_no_interaction(rho[i + 1], omega); 

		if (i < (N - 2))
		{
			A(i, i + 1) = -1 / pow(h, 2);
			A(i + 1, i) = -1 / pow(h, 2);
		}

	}

// Eigenvector matrix, R
	mat R = zeros<mat>(N - 1, N - 1);
	for (int i = 0; i<N - 1; i++)
	{
		R(i, i) = 1.0;
	}

// Finding eigenvalues and eigenvector using Armadillo/LAPACK
	start_2 = clock();
	vec eigval; mat eigvec;

	eig_sym(eigval, eigvec, A);
	finish_2 = clock();

	int kmax;int lmax;
	int max_iteration = 5000000;
	int iteration = 0;
	double eps = pow(10, -9);

	start_1 = clock();

// Finding largest non-diagonal element
	largest_element(A, &kmax, &lmax);
// Finding eigenvalues
	while (abs(A(kmax, lmax)) > eps && iteration < max_iteration)
	{
		iteration += 1;
		Jac_Rot(A, R, kmax, lmax, N);
		largest_element(A, &kmax, &lmax);
	}
	finish_1 = clock();

	cout << "Iterations:" << iteration << endl; // Printing number of iterations

	vec lambda = diagvec(A); // Saving eigenvalues

	cout << "CPU time, program: "  << double(finish_1 - start_1) / double(CLOCKS_PER_SEC) << "s" << endl;
	cout << "CPU time, armadillo: " << double(finish_2 - start_2) / double(CLOCKS_PER_SEC) << "s" << endl;
	
	ortho_test(R, eps); // Testing orthogonality
	test_largest_element(2, 5); // Testing the largest_element function

// Boundary conditions
	mat u = zeros<mat>(N + 1, N - 1);
	for (int i = 0; i < N - 1; i++)
	{
		for (int j = 0; j < N - 1; j++)
		{
			u(i + 1, j) = R(i, j);
		}
	}



// Algorithm for finding the three lowest eigenvalues
	int k1 = 0; int k2 = 0; int k3 = 0;
	double lambda1 = 500.0; double lambda2 = 500.0; double lambda3 = 500.0; 

	for (int k = 0; k < N - 1; k++)
	{
		if (lambda[k] < lambda1)
		{
			lambda3 = lambda2;
			lambda2 = lambda1;
			lambda1 = lambda[k];
			k3 = k2; k2 = k1; k1 = k;
		}
		if (lambda[k]  < lambda2 && lambda[k] > lambda1)
		{
			lambda3 = lambda2;
			lambda2 = lambda[k];
			k3 = k2; k2 = k;
		}
		if (lambda[k]  < lambda3 && lambda[k] > lambda2)
		{
			lambda3 = lambda[k];
			k3 = k;
		}

	}


// Printing eigenvalues
		cout << "Position: "<< k1 << "  " << k2 << "  " << k3 << endl;
		cout << "Eigenvalues: " << lambda[k1] << "  " << lambda[k2] << "  " << lambda[k3] << endl;

// Writing eigenvalues to file
	ofstream eigenvalues;
	eigenvalues.open("../eigenvalues_2e_no_int.txt");
	for (int i = 0; i < lambda.n_rows; i++)
	{
		eigenvalues << lambda[i] << "  ";
	}
	eigenvalues.close();

// Eigenvectors
	vec eigenvec1 = u.col(k1); vec eigenvec2 = u.col(k2); vec eigenvec3 = u.col(k3);

// Writing eigenvectors to file
	ofstream eigenvectors;
	eigenvectors.open("../eigenvectors_2e_no_int.txt");
    for(int i=0; i < N+1; i++){
		eigenvectors << eigenvec1[i] << " " << eigenvec2[i] << " " << eigenvec3[i] << endl;
	}
	eigenvectors.close();

	getchar();
}