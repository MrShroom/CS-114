#include "MainForPart1Task3.h"

namespace Task3 {

    int main_2(int argc, char *argv[])
    {
        MonteCarloIntegrationSolver *mySolver = new MonteCarloIntegrationSolver();
        std::fstream myfile;
        myfile.open("outputPartOne.txt", std::fstream::out | std::fstream::app);
        
        //x,y,z,theta, phi
        double a[] = { -0.5,-0.5, 0.0, 0.0, 0.0 };
        double b[] = { 0.5, 0.5, 0.0, 1.0, 1.0 };

        mySolver->setFuntion(&taskThreeFuntion);
        mySolver->setDensityFunctionAsUni(&taskThreeDensityFuntion, a, b, 5);

        myfile << "\nTask Three:\n";
        outTitleBar(myfile);
        SolverResults *results;
        int NumberOfSamples = 2000000;
        double current_s_bar;

        do
        {
            mySolver->setNumberOfSamples(NumberOfSamples);
            results = mySolver->solveArr(5);
            current_s_bar = sovleForSBar(results);
            myfile << "| " << std::setw(8) << std::left << NumberOfSamples;
            myfile << "| " << std::setw(8) << std::left << std::setprecision(3) << std::setiosflags(std::ios::fixed) << results->I_bar;
            myfile << "| " << std::setw(8) << std::left << std::setprecision(3) << std::setiosflags(std::ios::fixed) << current_s_bar;
            myfile << std::setw(1) << "|";
            myfile << std::endl;

            if (current_s_bar > 0.15)
                NumberOfSamples += 1500000;
            else if (current_s_bar > 0.125)
                NumberOfSamples += 1000000;
            else  if (current_s_bar > 0.1125)
                NumberOfSamples += 500000;
            else if(current_s_bar > 0.105)
                NumberOfSamples += 100000;
            else if(current_s_bar > 0.101)
                NumberOfSamples += 50000;
            else if (current_s_bar > 0.1001)
                NumberOfSamples += 10000;
            else 
                NumberOfSamples += 1000;

            delete(results);
        } while (current_s_bar > 0.1);

        return 0;
    }

    double chi(double * ray_coord, double * omega, double * light_center, double radius_of_light)
    {
        //convert into cartisain vector
        double v[3];
        v[0] = sin(omega[0]) * cos(omega[1]);
        v[1] = sin(omega[0]) * sin(omega[1]);
        v[2] = cos(omega[0]);

        //translate point to center at orgin
        double w[3];
        w[0] = light_center[0] - ray_coord[0];
        w[1] = light_center[1] - ray_coord[1];
        w[2] = light_center[2] - ray_coord[2];

        //find closest point on vector to shpere center
        double dist_from_ray_cood =
            (w[0] * v[0] + w[1] * v[1] + w[2] * v[2]) /
            (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

        double point_on_ray[3];
        point_on_ray[0] = ray_coord[0] + dist_from_ray_cood * v[0];
        point_on_ray[1] = ray_coord[1] + dist_from_ray_cood * v[1];
        point_on_ray[2] = ray_coord[2] + dist_from_ray_cood * v[2];

        //get squared ditance to light center
        double d =
            pow(point_on_ray[0] - light_center[0], 2) +
            pow(point_on_ray[1] - light_center[1], 2) +
            pow(point_on_ray[2] - light_center[2], 2);

        if (d > radius_of_light*radius_of_light)
            return 0.0;
        return 1.0;
    }

    double taskThreeFuntion(double * x, int size)
    {

        double light_center[3] = { 1.0,1.0,5.0 };
        // x[3] = theta, x[4] phi
        double omega[2];

        omega[0] = acos(x[3]);
        omega[1] = 2.0 * PI * x[4];
        return 100.0 * chi(x, omega, light_center, 1.0) * cos(omega[0]);
    }

    double taskThreeDensityFuntion(double * x, int size)
    {
        return 1.0 / (2.0*PI);
    }

    double sovleForSBar(SolverResults * results)
    {
        double s_bar = 0;
        for (size_t i = 0; i < results->I.size(); i++)
        {
            s_bar += pow((results->I[i] - results->I_bar), 2.0);
        }
        s_bar /= (results->I.size() - 1.0);
        s_bar = sqrt(s_bar);
        s_bar *= 2;
        s_bar /= sqrt(results->I.size());
        return s_bar;
    }

    void outTitleBar(std::fstream & myfile)
    {
        myfile << std::setw(10) << std::left << "| # Samples:";
        myfile << std::setw(10) << std::left << "| I_bar";
        myfile << std::setw(10) << std::left << "| s_bar";
        myfile << std::setw(1) << std::left << "|";
        myfile << std::endl;
    }
};