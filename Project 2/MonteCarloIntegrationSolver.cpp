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

		results->I[i] = m_funtion_solve(x) / m_funtion_density(x);
		sum_of_I += results->I[i];
	}
	results->I_bar = sum_of_I / NUMBER_OF_SAMPLES;
	return results;
}

SolverResults_t * MonteCarloIntegrationSolver::solveArr(int arraySize)
{
    SolverResults_t* results = new SolverResults_t();
    double *x = new double[arraySize];
    double sum_of_I = 0.0;

    for (int i = 0; i < NUMBER_OF_SAMPLES; i++)
    {
        for (int j = 0; j < arraySize; j++)
        {
            switch (m_curreunt_distribution)
            {
            case UNIFORM:
                x[j] = m_uniform_distribution[j](*m_generator);
                break;
            case EXPONENTENIAL:
                x[j] = m_exponential_distribution[j](*m_generator);
                break;
            default:
                x[j] = 0;
                break;
            }
        }
        results->I[i] = m_funtion_solve_arr(x, arraySize) / m_funtion_density_arr(x, arraySize);
        sum_of_I += results->I[i];
    }
    results->I_bar = sum_of_I / NUMBER_OF_SAMPLES;
    return results;
}


void MonteCarloIntegrationSolver::setFuntion(function_ptr funtion_solve)
{
	m_funtion_solve = funtion_solve;
}

void MonteCarloIntegrationSolver::setFuntion(function_ptr_array function_solve_array)
{
    m_funtion_solve_arr = function_solve_array; 
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

void MonteCarloIntegrationSolver::setDensityFunctionAsExp(function_ptr_array funtion_density, double lambda, int array_size)
{
    m_funtion_density_arr = funtion_density;
    m_curreunt_distribution = EXPONENTENIAL;
    if (m_exponential_distribution != 0)
    {
        delete(m_exponential_distribution);
    }
    m_exponential_distribution = new std::exponential_distribution<double>[array_size];
    for (int i = 0; i < array_size; i++)
    {
        m_exponential_distribution[i] = std::exponential_distribution<double>(lambda);
    }
}

void MonteCarloIntegrationSolver::setDensityFunctionAsUni(function_ptr_array funtion_density, double a[], double b[] ,int array_size)
{
    m_funtion_density_arr = funtion_density;
    m_curreunt_distribution = UNIFORM;
    if (m_uniform_distribution != 0)
    {
        delete(m_uniform_distribution);
    }
    m_uniform_distribution = new std::uniform_real_distribution<double>[array_size];
    for (int i = 0; i < array_size; i++)
    {
        m_uniform_distribution[i] = std::uniform_real_distribution<double>(a[i], b[i]);
    }
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

