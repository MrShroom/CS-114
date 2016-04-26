#pragma once
#include<random>

#define NUMBER_OF_SAMPLES 100000

typedef double(*function_ptr)(double);

enum Distro
{
	UNIFORM,
	EXPONENTENIAL
};

typedef struct SolverResults {
	double I[NUMBER_OF_SAMPLES];
	double I_bar;
} SolverResults_t ;

class MonteCarloIntegrationSolver {
public:
	MonteCarloIntegrationSolver();
	SolverResults_t* solve();
	void setFuntion(function_ptr funtion_slove);
	void setDensityFunctionAsExp(function_ptr funtion_density, double lambda);
	void setDensityFunctionAsUni(function_ptr funtion_density, double a, double b);
	~MonteCarloIntegrationSolver();
private:
	function_ptr m_funtion_slove;
	function_ptr m_funtion_density;
	std::random_device rd;
	std::mt19937 *m_generator = new std::mt19937(rd());
	std::uniform_real_distribution<double> *m_uniform_distribution = 0;
	std::exponential_distribution<double> *m_exponential_distribution = 0;
	Distro m_curreunt_distribution;
};