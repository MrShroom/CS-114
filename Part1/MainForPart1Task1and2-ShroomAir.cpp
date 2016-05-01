#include "MainForPArt1Task1and2.h"

namespace Task1and2 {

    extern double lamba_for_one_two = 0.1;

    int main_1(int argc, char *argv[])
    {
        MonteCarloIntegrationSolver *mySolver = new MonteCarloIntegrationSolver();

        //Start of Task 1-1
        std::ofstream myfile;
        myfile.open("outputPartOne.txt");
        myfile << "Task One-One:\n";
        outTitleBar(myfile);
        mySolver->setFuntion(&taskOneOneFuntion);
        mySolver->setDensityFunctionAsUni(&taskOneOneDensityFuntion, -2.0, 2.0);
        runSovlerTenTimes(mySolver, myfile);

        //Start of Task 1-2
        mySolver->setFuntion(&taskOneTwoFuntion);
        mySolver->setDensityFunctionAsUni(&taskOneTwoDensityFuntion, 0.0, 1.0);

        myfile << "\nTask One-Two:";
        for (int j = 0; j < 3; j++)
        {
            myfile << "\nlambda = " << lamba_for_one_two << std::endl;
            outTitleBar(myfile);
            runSovlerTenTimes(mySolver, myfile);
            lamba_for_one_two *= 10;
        }

        //Start of Task 2-1

        double a[] = { 0.0,0.0 };
        double b[] = { 1.0,1.0 };

        mySolver->setFuntion(&taskTwoOneFuntion);
        mySolver->setDensityFunctionAsUni(&taskTwoOneDensityFuntion, a, b, 2);

        myfile << "\nTask Two-One:\n";

        outTitleBar(myfile);
        runMultiDSovlerTenTimes(mySolver, myfile);


        //Start of Task 2-2

        mySolver->setFuntion(&taskTwoTwoFuntion);
        mySolver->setDensityFunctionAsUni(&taskTwoTwoDensityFuntion, a, b, 2);

        myfile << "\nTask Two-Two:\n";

        outTitleBar(myfile);
        runMultiDSovlerTenTimes(mySolver, myfile);

        myfile.close();

        delete(mySolver);
        return 0;
    }

    double taskOneOneFuntion(double x)
    {
        return exp(-pow(x, 2) / 2.0);
    }

    double taskOneOneDensityFuntion(double x)
    {
        if (x >= -2 && x < 2)
            return 0.25;
        return 0;
    }

    double taskOneTwoFuntion(double x)
    {
        double new_x = 1.0 - (log(x) / lamba_for_one_two);
        return taskOneOneFuntion(new_x);
    }

    double taskOneTwoDensityFuntion(double x)
    {
        double new_x = 1.0 - (log(x) / lamba_for_one_two);
        return lamba_for_one_two * exp(-lamba_for_one_two * (new_x - 1.0));
    }

    double taskTwoOneFuntion(double * x, int size)
    {
        return sin(x[0] * PI / 2.0);
    }

    double taskTwoOneDensityFuntion(double * x, int size)
    {
        return 1 / pow(PI, 2);
    }

    double taskTwoTwoFuntion(double * x, int size)
    {
        double new_x[2];
        new_x[0] = acos(x[0]);
        new_x[1] = 2.0 * PI * x[1];
        double out = cos(new_x[0]) * sin(new_x[0]) * cos(new_x[1]);
        return pow(out, 2);
    }

    double taskTwoTwoDensityFuntion(double * x, int size)
    {
        return 1.0 / (PI * PI);
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

    void outTitleBar(std::ofstream & myfile)
    {
        myfile << std::setw(10) << std::left << "| round:";
        myfile << std::setw(10) << std::left << "| I_bar";
        myfile << std::setw(10) << std::left << "| s_bar";
        myfile << std::setw(1) << std::left << "|";
        myfile << std::endl;
    }

    void runSovlerTenTimes(MonteCarloIntegrationSolver *solver, std::ofstream &myfile)
    {
        SolverResults *results;
        for (int i = 0; i < 10; i++)
        {
            results = solver->solve();
            myfile << "| " << std::setw(8) << std::left << i + 1;
            myfile << "| " << std::setw(8) << std::left << std::setprecision(3) << std::setiosflags(std::ios::fixed) << results->I_bar;
            myfile << "| " << std::setw(8) << std::left << std::setprecision(3) << std::setiosflags(std::ios::fixed) << sovleForSBar(results);
            myfile << std::setw(1) << "|";
            myfile << std::endl;
            delete(results);
        }
    }

    void runMultiDSovlerTenTimes(MonteCarloIntegrationSolver * sovler, std::ofstream & myfile)
    {
        SolverResults *results;
        for (int i = 0; i < 10; i++)
        {
            results = sovler->solveArr(2);
            myfile << "| " << std::setw(8) << std::left << i + 1;
            myfile << "| " << std::setw(8) << std::left << std::setprecision(3) << std::setiosflags(std::ios::fixed) << results->I_bar;
            myfile << "| " << std::setw(8) << std::left << std::setprecision(3) << std::setiosflags(std::ios::fixed) << sovleForSBar(results);
            myfile << std::setw(1) << "|";
            myfile << std::endl;
            delete(results);
        }
    }
}