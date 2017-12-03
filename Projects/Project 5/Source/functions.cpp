#include "functions.h"
using namespace std;

Setup::Setup(int N_cycles, int N_agents, double lambda, double alpha, double gamma, double money_initial, double bin_width) :
	m_cycles(N_cycles), m_money_initial(money_initial), m_agents(N_agents), m_lambda(lambda), m_alpha(alpha), m_gamma(gamma), m_bin_width(bin_width) {

	// Generate filenames based on input parameters
	stringstream filename;
	filename << setiosflags(ios::showpoint) << setprecision(3);
	filename << "MONEY_agents_" << N_agents << "_MCC_" << N_cycles << "_N_trans_" << m_transactions << "_lambda_" << lambda << "_alpha_" << alpha << "_gamma_" << gamma << ".dat";
	pdf_file = filename.str();

	stringstream filename2;
	filename2 << setiosflags(ios::showpoint) << setprecision(3);
	filename2 << "CYCLE_agents_" << N_agents << "_MCC_" << N_cycles << "_N_trans_" << m_transactions << "_lambda_" << lambda << "_alpha_" << alpha << "_gamma_" << gamma << ".dat";

	export_file.open(filename2.str().c_str());
	export_file << setiosflags(ios::showpoint) << setprecision(8);
	export_file << "cycle  " << "sigma  " << "expec  " << endl;

	m_bins = m_money_initial * m_agents / m_bin_width;
	bins_arr = new long[m_bins];

	for (int i = 0; i < m_bins; i++) {
		bins_arr[i] = 0;
	}
}

void Setup::start(mt19937_64& generate) {
	for (int i = 0; i < m_cycles; i++) {
		Calculation calculation(m_lambda, m_money_initial, m_agents, m_alpha, m_gamma);
		calculation.start(generate);
		bins_fix(calculation, i);
	}
}

void Setup::export_pdf() {
	ofstream histogram;
	histogram.open(pdf_file.c_str());
	histogram << setiosflags(ios::showpoint);

	// Export probability distribution
	for (int i = 0; i < m_bins; i++) {
		histogram << setprecision(9) << (m_bin_width*i) << " " << setprecision(8) << (bins_arr[i] / (double)m_cycles) << endl;
	}
	histogram.close();
}


void Setup::export_data_cycles(Calculation& calculation, int cycle) {
	// Export cycle data (cycle number, variance and expectation value)
	double sigma = 0;  double expec = 0;

	for (int i = 0; i < m_bins; i++) {
		double m = m_bin_width*(i + 0.5);
		double m2 = pow(m, 2);
		double pdf = bins_arr[i] / (double)cycle / (double)m_agents;

		sigma += pdf*m2;
		expec += pdf*m;
	}

	sigma -= pow(expec, 2);
	export_file << cycle << " " << sigma << " " << expec << endl;
}


void Setup::bins_fix(Calculation& calculation, int cycle) {
	double* money = calculation.get_wealth();
	for (int i = 0; i < m_agents; i++) {
		int bin = (int)floor(money[i] / m_bin_width);
		bins_arr[bin]++;
	}

	if (cycle > 1) {
		export_data_cycles(calculation, cycle);
	}
}


Calculation::Calculation(double lambda, double start_money, int agents, double alpha, double gamma) :
	m_agents(agents), m_lambda(lambda), m_alpha(alpha), m_gamma(gamma) {

	m_wealth = new double[m_agents];
	m_c = (int**)matrix(m_agents, m_agents, sizeof(int));
	m_interactions = 0;
	m_interactions_max = 0;
	for (int i = 0; i < m_agents; i++) {
		m_wealth[i] = start_money;
		for (int j = 0; j < m_agents; j++) {
			m_c[i][j] = 0;
		}
	}
}


void Calculation::start(mt19937_64& generate) {
	// Picking random numbers from uniform distributions from std
	uniform_real_distribution<double> choose_amount(0.0, 1.0);
	uniform_real_distribution<double> choose_interaction(0.0, 1.0);
	uniform_int_distribution<int> choose_agent(0, m_agents - 1); 
	

	// Main algorithm
	for (int i = 0; i < m_transactions; i++) {
		int a = choose_agent(generate); int b = choose_agent(generate);
		double prob;
		double eps = choose_amount(generate);
		double m_diff = fabs(m_wealth[a] - m_wealth[b]);
		
		if ((fabs(m_alpha) < 1e-8) || m_diff < 1){
			prob = 1;
		}
		else {
			prob = pow(m_diff, -m_alpha) * pow( (m_c[a][b] + 1), m_gamma);
		}

		double interaction = choose_interaction(generate);

		if (prob >= interaction) {
			double delta_m = (1 - m_lambda)*(eps*m_wealth[b] - (1 - eps)*m_wealth[a]);
			m_wealth[a] += delta_m;
			m_wealth[b] -= delta_m;

			m_c[a][b] = m_c[a][b] + 1;
			m_c[b][a] = m_c[a][b];

			if (m_interactions_max < m_c[a][b]) {
				m_interactions_max = m_c[a][b];
				m_interactions++;
			}
		}
	}
}








