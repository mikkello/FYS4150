#include "functions.h"
#include <ctime>

using namespace std;

int main() {

	// Parameters
	int N_cycles = 1e3; // MC cycles
	int N_agents = 1000; // Number of agents
	double money_initial = 1.0; // Initial money
	double bin_width = 0.05; // Bin width for histogram

	double lambda = 0.5; // Savings rate
	double alpha = 2.0; // Nearest neighbour parameter p_ij \propto |m_i - m_j|^(-alpha)
	double gamma = 4.0; // Memory parameter p_ij \propto (c_ij + 1)^(gamma)
	//double min_gamma = 0.0;
	//double max_gamma = 4.0;
	//double step_gamma = 1.0;

	// Random number generator (RNG) seed
	random_device seed;
	mt19937_64 generate(seed());
	clock_t t_start, t_finish;

	// Start calculation
	//for (double gamma = min_gamma; gamma <= max_gamma; gamma += step_gamma) {
	t_start = clock();
	Setup setup(N_cycles, N_agents, lambda, alpha, gamma, money_initial, bin_width);
	setup.start(generate);

	// Save distribution to file
	setup.export_pdf();

	// Timing
	t_finish = clock();
	double time = (t_finish - t_start);
	cout << "Calculations done." << endl;
	cout << "CPU time: " << (time / CLOCKS_PER_SEC) << "s" << endl;
	//}
}