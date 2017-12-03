#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <random>
#include <cmath>
#include "lib.h"


class Calculation {
private:
	int m_agents;
	int m_interactions;
	int m_interactions_max;

	int m_transactions = 1e7;

	double m_lambda;
	double m_alpha;
	double m_gamma;
	double* m_wealth;

	int** m_c;
	


public:
	Calculation(double lambda, double start_money, int N_agents, double alpha, double gamma);

	double* get_wealth() { return m_wealth; }
	int interactions() { return m_interactions; }
	int interactions_max() { return m_interactions_max; }
	void start(std::mt19937_64& generate);
	

};

class Setup {
private:
	int m_cycles;
	double m_money_initial;
	int m_agents;
	double m_bin_width;

	double m_lambda;
	double m_alpha;
	double m_gamma;
	

	int m_bins;
	long* bins_arr;

	int m_transactions = 1e7;

	std::string pdf_file;
	std::ofstream export_file;

public:

	Setup(int N_cycles, int N_agents, double lambda, double alpha, double gamma, double money_initial , double bin_width);

	void start(std::mt19937_64& generate);
	void bins_fix(Calculation& calculation, int cycle);

	void export_data_cycles(Calculation& calculation, int cycle);
	void export_pdf();
};