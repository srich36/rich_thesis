#include <iostream>
#include <string>
#include <math.h>

using namespace std;

double testParticle[] = {-0.856580908970521,-0.928947964866100,-0.532194537192831,-0.410209473732025,-0.0110478865577260,-0.654732028329036,0.937910086713157,-0.964208556080876,1.91880746165091,4.43391078133062,0.581192538349102 };
//double testParticle[] = {    -1.0000    0.1971    1.0000    0.9083   -0.3037   -1.0000    0.4753   -1.0000    0.7950    2.6223    0.4526  };

/*
************Configuration params************
*/
bool evalTestParticle = false;
bool debug = false;
bool debugPSO = false;

/*
************Rehydration params************
*/
bool rehydrate = false;
double averageJThresholdPercent = .1;
int numIterBeforeHydration = 100;
int numPrevIterationsAvgJ = 50;
double rehydratePercentage = .33;

/*
************Display params************
*/
bool displayCostFunctionValue = false;
bool displayAverageCostValuePerIteration = false;
bool displayGlobalBestPerIteration = true;
bool displayIterationNum = false;
bool displayExecutionNum = true;
bool displayPsoIterationNumber = true;
bool displayGlobalBestPerPSOIteration = true;

/*
************PSO Penalty params************
*/
double badResultPenalty = 100000;
int penaltyCoefficient = 100;
double penaltyValueCutoff = .001;


struct fuelParams{
    double c;
    double n0;
    fuelParams(double c=.5, double n0 = .2){
        this->c = c;
        this->n0 = n0;
    }
};

/* Eventually this should read all the params from a config file
   so we don't have to rebuild it everytime we modify a param
*/

int main(){
    /*PSO Configuration params*/
    int Beta = 2;

    /* Initial condition set */
    double vrInitial = 0; double vrTerminal = 0; double xiInital = 0;
    double R1 = 1; double R2 = Beta*R1; double ub = 1;
    double rInitial = R1; double rFinal = R2;
    double vThetaInital = sqrt(ub/R1); double vThetaFinal = sqrt(ub/R2);

    double cFuel = .5;
    double n0Fuel = .2;
    fuelParams fuel = fuelParams(cFuel, n0Fuel);

    cout << fuel.c << endl;

    return 0;
}