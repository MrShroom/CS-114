#pragma once
#ifndef MAIN_PART_1_AND_2_H
#define MAIN_PART_1_AND_2_H

#include "MonteCarloIntegrationSolver.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif // !PI

namespace Task1and2 {
    extern double lamba_for_one_two;

    double taskOneOneFuntion(double x);
    double taskOneOneDensityFuntion(double x);

    double taskOneTwoFuntion(double x);
    double taskOneTwoDensityFuntion(double x);

    double taskTwoOneFuntion(double *x, int size);
    double taskTwoOneDensityFuntion(double *x, int size);

    double taskTwoTwoFuntion(double *x, int size);
    double taskTwoTwoDensityFuntion(double *x, int size);

    double sovleForSBar(SolverResults * results);

    void outTitleBar(std::ofstream &myfile);

    void runSovlerTenTimes(MonteCarloIntegrationSolver *sovler, std::ofstream &myfile);
    void runMultiDSovlerTenTimes(MonteCarloIntegrationSolver *sovler, std::ofstream &myfile);

    int main_1(int argc, char *argv[]);
};

#endif // !MAIN_PART_1_AND_2_H