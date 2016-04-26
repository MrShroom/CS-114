#include "MonteCarloIntegrationSolver.h"


MonteCarloIntegrationSolver::MonteCarloIntegrationSolver()
{
}

SolverResults_t* MonteCarloIntegrationSolver::solve()
{
	SolverResults_t* results = new SolverResults_t();
	double x = 0.0;
	double sum_of_I = 0.0;
	
	for (int i = 0; i < NUMBER_OF_SAMPLES; i++)
	{
		switch (m_curreunt_distribution)
		{
		case UNIFORM:
			x= (*m_uniform_distribution)(*m_generator);
			break;
		case EXPONENTENIAL:
			x = (*m_exponential_distribution)(*m_generator);
			break;
		default:
			x = 0;
			break;
		}

		results->I[i] = m_funtion_slove(x) / m_funtion_density(x);
		sum_of_I += results->I[i];
	}
	results->I_bar = sum_of_I / NUMBER_OF_SAMPLES;
	return results;
}


void MonteCarloIntegrationSolver::setFuntion(function_ptr funtion_slove)
{
	m_funtion_slove = funtion_slove;
}

void MonteCarloIntegrationSolver::setDensityFunctionAsExp(function_ptr funtion_density, double lambda)
{
	m_funtion_density = funtion_density;
	m_curreunt_distribution = EXPONENTENIAL;
	if (m_exponential_distribution != 0)
	{
		delete(m_exponential_distribution);
	}
	m_exponential_distribution = new std::exponential_distribution<double>(lambda);
}

void MonteCarloIntegrationSolver::setDensityFunctionAsUni(function_ptr funtion_density, double a, double b)
{
	m_funtion_density = funtion_density;
	m_curreunt_distribution = UNIFORM;
	if (m_uniform_distribution != 0)
	{
		delete(m_uniform_distribution);
	}
	m_uniform_distribution = new std::uniform_real_distribution<double>(a, b);
}

MonteCarloIntegrationSolver::~MonteCarloIntegrationSolver()
{
	if (m_generator != 0)
		delete(m_generator);
	if (m_uniform_distribution != 0)
		delete(m_uniform_distribution);
	if (m_exponential_distribution != 0)
		delete(m_exponential_distribution);
}

