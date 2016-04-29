#pragma once
#ifndef MAIN_FOR_PART_1_TASK_3_H
#define MAIN_FOR_PART_1_TASK_3_H

#include "MonteCarloIntegrationSolver.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

namespace Task3 {

    double chi(double *ray_coord, double *omega, double *light_center, double radius_of_light);

    double taskThreeFuntion(double *x, int size);
    double taskThreeDensityFuntion(double *x, int size);

    void outTitleBar(std::fstream &myfile);
    double sovleForSBar(SolverResults * results);

    int main_2(int argc, char *argv[]);
}
#endif // !MAIN_FOR_PART_1_TASK_3_H