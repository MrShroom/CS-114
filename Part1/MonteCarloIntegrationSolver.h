#pragma once
#include <random>
#include <vector>

typedef double(*function_ptr)(double);
typedef double(*function_ptr_array)(double* , int );

enum Distro
{
	UNIFORM,
	EXPONENTENIAL
};

typedef struct SolverResults {
	std::vector<double> I;
	double I_bar;
} SolverResults_t ;

class MonteCarloIntegrationSolver {
public:
	MonteCarloIntegrationSolver();
	SolverResults_t* solve();
    SolverResults_t* solveArr(int arraySize);

	void setFuntion(function_ptr funtion_slove);
    void setFuntion(function_ptr_array function_solve_array);
    void setNumberOfSamples(int number_of_samples);

	void setDensityFunctionAsExp(function_ptr funtion_density, double lambda);
	void setDensityFunctionAsUni(function_ptr funtion_density, double a, double b);

    void setDensityFunctionAsExp(function_ptr_array funtion_density, double lambda, int arraySize);
    void setDensityFunctionAsUni(function_ptr_array funtion_density, double a[], double b[], int array_size);

	~MonteCarloIntegrationSolver();
private:
	function_ptr m_funtion_solve;
    function_ptr_array m_funtion_solve_arr;
	function_ptr m_funtion_density;
    function_ptr_array m_funtion_density_arr;
    
    int  m_number_of_samples = 100000;
	std::random_device rd;
	std::mt19937 *m_generator = new std::mt19937(rd());
	std::uniform_real_distribution<double> *m_uniform_distribution = 0;
	std::exponential_distribution<double> *m_exponential_distribution = 0;
	Distro m_curreunt_distribution;
};