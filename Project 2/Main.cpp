#include "MonteCarloIntegrationSolver.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

double lamba_for_one_two = 0.1;

double taskOneOneFuntion(double x);
double taskOneOneDensityFuntion(double x);

double taskOneTwoFuntion(double x);
double taskOneTwoDensityFuntion(double x);

double sovleForSBar(SolverResults * results);

int main()
{
	MonteCarloIntegrationSolver *mySolver = new MonteCarloIntegrationSolver();
	
	//Start of Task 1-1
	std::ofstream myfile;
	myfile.open("outputPartOne.txt");
	myfile << "Task One-One:\n";
	myfile << std::setw(10) << std::left << "| round:";
	myfile << std::setw(10) << std::left << "| I_bar";
	myfile << std::setw(10) << std::left << "| s_bar";
	myfile << std::setw(1) << std::left << "|";
	myfile << std::endl;

	mySolver->setFuntion(&taskOneOneFuntion);
	mySolver->setDensityFunctionAsUni(&taskOneOneDensityFuntion, -2.0, 2.0);

	SolverResults *results;
	for (int i = 0; i < 10; i++)
	{
		results = mySolver->solve();
		
		myfile << "| " << std::setw(8) << std::left  << i + 1;;
		myfile << "| "<< std::setw(8) << std::left << std::setprecision(3)  << std::setiosflags(std::ios::fixed)
			<<  results->I_bar;
		myfile << "| " << std::setw(8) << std::left << std::setprecision(3) << std::setiosflags(std::ios::fixed)
			<< sovleForSBar(results);
		myfile << std::setw(1) << "|";
		myfile << std::endl;
		delete(results);
	}

	//Start of Task 1-2
	mySolver->setFuntion(&taskOneTwoFuntion);
	mySolver->setDensityFunctionAsUni(&taskOneTwoDensityFuntion, 0.0, 1.0);

	myfile << "Task One-Two:\n";
	for (int j = 0; j < 3; j++)
	{
		myfile << "\nlambda = " << lamba_for_one_two << std::endl;
		myfile << std::setw(10) << std::left << "| round:";
		myfile << std::setw(10) << std::left << "| I_bar";
		myfile << std::setw(10) << std::left << "| s_bar";
		myfile << std::setw(1) << std::left << "|";
		myfile << std::endl;

		for (int i = 0; i < 10; i++)
		{
			results = mySolver->solve();

			myfile << "| " << std::setw(8) << std::left << i + 1;;
			myfile << "| " << std::setw(8) << std::left << std::setprecision(3) << std::setiosflags(std::ios::fixed)
				<< results->I_bar;
			myfile << "| " << std::setw(8) << std::left << std::setprecision(3) << std::setiosflags(std::ios::fixed)
				<< sovleForSBar(results);
			myfile << std::setw(1) << "|";
			myfile << std::endl;
			delete(results);
		}
		lamba_for_one_two *= 10;
	}
	
	//Start of Task 2-1

	myfile.close();
	
	delete(mySolver);
	return 0;
}

double taskOneOneFuntion(double x)
{
	return exp(-pow(x, 2) / 2);
}

double taskOneOneDensityFuntion(double x)
{
	if (x >= -2 && x < 2)
		return 0.25;
	return 0;
}

double taskOneTwoFuntion(double x)
{
	double new_x = 1 - (log(x) / lamba_for_one_two);
	return taskOneOneFuntion(new_x);
}

double taskOneTwoDensityFuntion(double x)
{
	double new_x = 1 - (log(x) / lamba_for_one_two);
	return lamba_for_one_two * exp(- lamba_for_one_two * (new_x -1.0));
}

double sovleForSBar(SolverResults * results)
{
	double s_bar = 0;
	for (int i = 0; i < NUMBER_OF_SAMPLES; i++)
	{
		s_bar += pow((results->I[i] - results->I_bar), 2.0);
	}
	s_bar /= (NUMBER_OF_SAMPLES - 1.0);
	s_bar = sqrt(s_bar);
	s_bar *= 2;
	s_bar /= sqrt(NUMBER_OF_SAMPLES);
	return s_bar;
}
