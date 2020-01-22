#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <vector>
#include <stdlib.h>
#include <time.h>

using namespace std;


/*
************Configuration params************
*/
const int numExecutions=1, numPsoIterations=1, numParticles = 100;
const int numUnknowns = 11, numIterations = 500;
/*
************Configuration params************
*/

struct fuelParams{
    double c;
    double n0;
    fuelParams(double c=.5, double n0 = .2){
        this->c = c;
        this->n0 = n0;
    }
};


class particleBounds{
    private:
        double t1, t2, xi, v, de;
        double bounds[numUnknowns];
    public:
        particleBounds(double t1, double de, double t2, double xi, double v){
            this->t1 = t1;
            this->t2 = t2;
            this->xi = xi;
            this->v = v;
            this->de = de;
            bounds[0]=xi; bounds[1] = xi; bounds[2]=xi; bounds[3] = xi;
            bounds[4]=v; bounds[5] = v; bounds[6]=v; bounds[7] = v;
            bounds[8] = t1; bounds[9] = de; bounds[10] = t2;
        }

        double* getBounds() { return bounds; }

};

struct thrustArcIc{
    double vrInitial;
    double vThetaInital;
    double rInitial;
    double xiInital;
    thrustArcIc(double vr, double vt, double r, double xi){
        vrInitial = vr;
        vThetaInital = vt;
        rInitial = r;
        xiInital = xi;
    }
};

void allocateSwarm(double *swarm, int numParticles, int numUnknowns, particleBounds, particleBounds );
void printSwarm(double *swarm, int, int);
void printParticleBests(double *, int);
int indexConversion(int, int, int);
double randNumberPSO(double, double);

double testParticle[] = {-0.856580908970521,-0.928947964866100,-0.532194537192831,-0.410209473732025,-0.0110478865577260,-0.654732028329036,0.937910086713157,-0.964208556080876,1.91880746165091,4.43391078133062,0.581192538349102 };
//double testParticle[] = {    -1.0000    0.1971    1.0000    0.9083   -0.3037   -1.0000    0.4753   -1.0000    0.7950    2.6223    0.4526  };



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

/*
************Parameter Bounds Params************
*/
const double LBt1=0.001, LBdeltaE = 0, LBdeltat2=0.001, LBxi = -1, LBv=-1;
const double UBt1=3, UBdeltaE = 2*M_PI, UBdeltat2=3, UBxi = 1, UBv=1;

/* Eventually this should read all the params from a config file
   so we don't have to rebuild it everytime we modify a param
*/

int main(){
    /*PSO Configuration params*/
    srand(time(0));
    cout << rand() / RAND_MAX << endl;
    int Beta = 2;

    /* Initial condition set */
    double vrInitial = 0; double vrTerminal = 0; double xiInital = 0;
    double R1 = 1; double R2 = Beta*R1; double ub = 1;
    double rInitial = R1; double rFinal = R2;
    double vThetaInital = sqrt(ub/R1); double vThetaFinal = sqrt(ub/R2);

    double cFuel = .5;
    double n0Fuel = .2;

    fuelParams fuel = fuelParams(cFuel, n0Fuel);
    thrustArcIc tarc1 = thrustArcIc(vrInitial, vThetaInital, rInitial, xiInital);
    particleBounds particleLB = particleBounds(LBt1, LBdeltaE, LBdeltat2, LBxi, LBv);
    particleBounds particleUB = particleBounds(UBt1, UBdeltaE, UBdeltat2, UBxi, UBv);

    for(int executionNum=0; executionNum<numExecutions; executionNum++){
        /*
        ************Main algorithm goes here************
        */
        double globalBestValue = numeric_limits<double>::max();
        vector<double> globalBestValuePerIteration;
        double globalBestParticle[numUnknowns];

        /*
            Allocate particles
        */
        double* swarm = new double[numParticles*numUnknowns];
        allocateSwarm(swarm, numParticles, numUnknowns, particleLB, particleUB);
        //Allocate Personal Particle Values
        double *particlePersonalBestValues = new double[numParticles];
        for(int i = 0; i < numParticles; i++){
            particlePersonalBestValues[i] = numeric_limits<double>::max();
        }

        double **particlePersonalBestParticles = new double*[numParticles];
        for(int i = 0; i < numParticles; i++){
            particlePersonalBestParticles[i] = new double[numUnknowns];
        }


        //printParticleBests(particlePersonalBestValues, numParticles);
        printSwarm(swarm, numParticles, numUnknowns);

        for(int iterationNum = 0; iterationNum<numIterations; iterationNum++){
            /*
            ************Evaluate a single iteration within this loop************
            */
        }


    }

    return 0;
}


void printSwarm(double *swarm, int numParticles, int numUnknowns){
    string swarmPrintString = "\n***********SWARM***********\n";
    cout << swarmPrintString << endl;
    for(int i = 0;i<numParticles; i++){
        for(int j = 0; j<numUnknowns; j++){
            cout << swarm[indexConversion(i,j, numUnknowns)] << " ";
        }
        cout << endl;
    }
    cout << swarmPrintString << endl;
}

void printParticleBests(double* particleBests, int numParticles){
    string particleBestsString = "\n***************PARTICLE BESTS**********\n";
    cout << particleBestsString << endl;
    for(int i=0; i < numParticles; i++){
        cout << particleBests[i] <<" ";
    }
}

void allocateSwarm(double *swarm, int numParticles, int numUnknowns, particleBounds lowerBounds, particleBounds upperBounds){
    cout << "Allocating swarm" << endl;
    double* lb = lowerBounds.getBounds();
    double* ub = upperBounds.getBounds();
    for(int i = 0;i<numParticles; i++){
        for(int j = 0; j<numUnknowns; j++){
            //Rand should go between 0 and 1
            double rand = randNumberPSO(0.0, 1.0);
            double randBetweenBounds = (ub[j]-lb[j])*rand+lb[j];
            swarm[indexConversion(i,j, numUnknowns)] = randBetweenBounds;
        }
    }
}

double randNumberPSO(double min, double max){
    double r3 = min + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(max-min)));
    return r3;
}

//(3,4) will return 4*11+3 = 44
int indexConversion(int x, int y, int arrWidth=numUnknowns){
    return x * arrWidth + y;
}