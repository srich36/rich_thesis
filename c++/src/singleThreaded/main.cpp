#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <time.h>
#include <sys/time.h>

using namespace std;
using namespace boost::numeric::odeint;


typedef std::vector<double> state_type;
typedef boost::numeric::odeint::result_of::make_dense_output<
    runge_kutta_dopri5<state_type>>::type dense_stepper_type;

const double MIN_STEP_SIZE = 7.105427e-15;
bool debug = false, outputParticleParams = false, evalTestParticle = false, outputWarnings = false;
double testParticle[] = {-0.6356245844181813, 0.6749135584289404, -0.3782733531510278, -0.1791742193502524, 0.2948498719982566, -0.2533456929439651, -0.9806058557165864, -0.1882857942679574, 0.7291294746670177, 3.083940333165867, 0.4104569474969907 };

/*
************Configuration params************
*/
const int numExecutions = 30, numPsoIterations = 1;
const int numUnknowns = 11;
int numParticles = 110;
int numIterations = 500;
/*
************Configuration params************
*/

class MinStepSizeException : public std::exception
{
    double vr, vTheta, r, xi;

public:
    MinStepSizeException(double vr, double vTheta, double r, double xi)
    {
        this->vr = vr;
        this->vTheta = vTheta;
        this->r = r;
        this->xi = xi;
    }

    double getVr() { return vr; }
    double getvTheta() { return vTheta; }
    double getR() { return r; }
    double getXi() { return xi; }
};

struct fuelParams
{
    double c;
    double n0;
    fuelParams(double c = .5, double n0 = .2)
    {
        this->c = c;
        this->n0 = n0;
    }
};

class particleBounds
{
private:
    double t1, t2, xi, v, de;
    double bounds[numUnknowns];

public:
    particleBounds(double t1, double de, double t2, double xi, double v)
    {
        this->t1 = t1;
        this->t2 = t2;
        this->xi = xi;
        this->v = v;
        this->de = de;
        bounds[0] = xi;
        bounds[1] = xi;
        bounds[2] = xi;
        bounds[3] = xi;
        bounds[4] = v;
        bounds[5] = v;
        bounds[6] = v;
        bounds[7] = v;
        bounds[8] = t1;
        bounds[9] = de;
        bounds[10] = t2;
    }

    double *getBounds() { return bounds; }
    friend particleBounds operator-(const particleBounds &p1, const particleBounds &p2);
    friend ostream &operator<<(ostream &os, const particleBounds &p1);
};

particleBounds operator-(const particleBounds &p1, const particleBounds &p2)
{
    return particleBounds(p1.t1 - p2.t2, p1.de - p2.de, p1.t2 - p2.t2, p1.xi - p2.xi, p1.v - p2.v);
}

ostream &operator<<(ostream &os, const particleBounds &p1)
{
    os << p1.bounds[0] << ' ' << p1.bounds[1] << ' ' << p1.bounds[2] << ' ';
    os << p1.bounds[3] << ' ' << p1.bounds[4] << ' ' << p1.bounds[5] << ' ';
    os << p1.bounds[6] << ' ' << p1.bounds[7] << ' ' << p1.bounds[8] << ' ';
    os << p1.bounds[9] << ' ' << p1.bounds[10] << ' ' << endl;
    return os;
}

struct thrustArcIc
{
    double vrInitial;
    double vThetaInital;
    double rInitial;
    double xiInital;
    thrustArcIc(double vr, double vt, double r, double xi)
    {
        vrInitial = vr;
        vThetaInital = vt;
        rInitial = r;
        xiInital = xi;
    }
};

class timeStepObserver
{

private:
    double prevTime, currTime, timeStep;
    int numSteps;

public:
    timeStepObserver(double initialTimeStep)
    {
        this->prevTime = 0;
        this->currTime = 0;
        this->numSteps = 0;
        this->timeStep = initialTimeStep;
    }

    void operator()(const state_type &x, double t)
    {
        numSteps++;
        swap(prevTime, currTime);
        currTime = t;
        timeStep = currTime - prevTime;
        if (timeStep < MIN_STEP_SIZE && currTime != 0)
        {
            throw(MinStepSizeException(x[0], x[1], x[2], x[3]));
        }
    }
};

class thrustArc1EOM
{

public:
    double ub, xi0, xi1, xi2, xi3;
    fuelParams fuel;

    thrustArc1EOM(double ub, double xi0, double xi1, double xi2, double xi3, fuelParams fuel)
    {
        this->ub = ub;
        this->xi0 = xi0;
        this->xi1 = xi1;
        this->xi2 = xi2;
        this->xi3 = xi3;
        this->fuel = fuel;
    }

    void operator()(const state_type &x, state_type &dxdt, const double t)
    {

        //Delta is working correctly
        double delta = xi0 + xi1 * t + xi2 * pow(t, 2) + xi3 * pow(t, 3);

        double ratio = (fuel.c * fuel.n0) / (fuel.c - fuel.n0 * t);

        //This integrates correctly for a good particle

        dxdt[0] = -1 * (ub - x[2] * pow(x[1], 2)) / (pow(x[2], 2)) + ratio * sin(delta);
        dxdt[1] = -1 * x[0] * x[1] / x[2] + ratio * cos(delta);
        dxdt[2] = x[0];
        dxdt[3] = x[1] / x[2];
    }
};

class thrustArc2EOM
{

public:
    double ub, v0, v1, v2, v3, deltaT1, coastingTimeInterval, t2;
    fuelParams fuel;

    thrustArc2EOM(double ub, double v0, double v1, double v2, double v3, fuelParams fuel, double deltaT1, double coastingTimeInterval)
    {
        this->ub = ub;
        this->v0 = v0;
        this->v1 = v1;
        this->v2 = v2;
        this->v3 = v3;
        this->fuel = fuel;
        this->deltaT1 = deltaT1;
        this->coastingTimeInterval = coastingTimeInterval;
        this->t2 = deltaT1 + coastingTimeInterval;
    }

    void operator()(const state_type &x, state_type &dxdt, const double t)
    {
        //Delta is working correctly
        //cout << t << endl;
        //delta = v0+v1.*(t-this->t2)+v2.*(t-t2)^2+v3.*(t-t2)^3;
        double delta = v0 + v1 * (t - t2) + v2 * pow((t - t2), 2) + v3 * pow((t - t2), 3);

        double ratio = (fuel.c * fuel.n0) / (fuel.c - fuel.n0 * (t + deltaT1 - t2));

        //This integrates correctly for a good particle

        dxdt[0] = -1 * (ub - x[2] * pow(x[1], 2)) / (pow(x[2], 2)) + ratio * sin(delta);
        dxdt[1] = -1 * x[0] * x[1] / x[2] + ratio * cos(delta);
        dxdt[2] = x[0];
        dxdt[3] = x[1] / x[2];
    }
};

//abc
void allocateSwarm(double *swarm, int numParticles, int numUnknowns, particleBounds, particleBounds);
void printSwarm(double *swarm, int, int);
void printValues(double *, int, string, int);
int indexConversion(int, int, int);
double randNumberPSO(double, double);
void setParticle(double *, double *, int, int);
void setGlobalBestParticle(double *, double *, int, int);
double getWallTime();
double getCPUTime();
void writeToFile(std::ofstream &, double, double, double, double *, int, int, int);
void writeHeadersToFile(std::ofstream &);
bool fileExists(string);

//double testParticle[] = {    -1.0000    0.1971    1.0000    0.9083   -0.3037   -1.0000    0.4753   -1.0000    0.7950    2.6223    0.4526  };

bool debugPSO = false;

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
bool displayPSOParams = false;
bool displayCostFunctionAverages = false;

/*
************PSO Penalty params************
*/
double badResultPenalty = 100000;
int penaltyCoefficient = 100;
double penaltyValueCutoff = .001;

/*
************Parameter Bounds Params************
*/
const double LBt1 = 0.001, LBdeltaE = 0, LBdeltat2 = 0.001, LBxi = -1, LBv = -1;
const double UBt1 = 3, UBdeltaE = 2 * M_PI, UBdeltat2 = 3, UBxi = 1, UBv = 1;

/* Eventually this should read all the params from a config file
   so we don't have to rebuild it everytime we modify a param
*/

typedef std::vector<double> state_type;

//Tolerances
const double absoluteTolerance = 1.0e-12;
const double relTolerance = 1.0e-12;

//Initial step size
const double initialStepSize = .0001;

int main()
{
    /*PSO Configuration params*/

    std::ofstream resultsFile;
    string fileName = "STpNum"+to_string(numParticles)+"Inum"+to_string(numIterations)+".csv";

    bool fileExistsBool = fileExists(fileName);
    resultsFile.open(fileName, ios::app);
    if(!fileExistsBool){
        writeHeadersToFile(resultsFile);
    }

    resultsFile.precision(16);
    resultsFile.setf(ios::fixed);
    resultsFile.setf(ios::showpoint);


    srand(time(0));
    int Beta = 2;

    /* Initial condition set */
    double vrInitial = 0;
    double vrTerminal = 0;
    double xiInital = 0;
    double R1 = 1;
    double R2 = Beta * R1;
    double ub = 1;
    double rInitial = R1;
    double rFinal = R2;
    double vThetaInital = sqrt(ub / R1);
    double vThetaFinal = sqrt(ub / R2);

    double cFuel = .5;
    double n0Fuel = .2;

    fuelParams fuel = fuelParams(cFuel, n0Fuel);
    thrustArcIc tarc1 = thrustArcIc(vrInitial, vThetaInital, rInitial, xiInital);
    particleBounds particleLB = particleBounds(LBt1, LBdeltaE, LBdeltat2, LBxi, LBv);
    particleBounds particleUB = particleBounds(UBt1, UBdeltaE, UBdeltat2, UBxi, UBv);

    for (int executionNum = 0; executionNum < numExecutions; executionNum++)
    {
        /*
        ************Main algorithm goes here************
        */

        //asd
        double globalBestValue = numeric_limits<double>::max();

        double *globalBestValuePerIteration = new double[numIterations];
        double *globalBestParticle = new double[numUnknowns];

        /*Start timer here */
        double wall0 = getWallTime();
        double cpu0 = getCPUTime();

        /*
            Allocate particles
        */
        double *swarm;
        if (evalTestParticle)
        {
            swarm = testParticle;
            numParticles = 1;
            numIterations = 1;
        }
        else
        {
            swarm = new double[numParticles * numUnknowns];
            allocateSwarm(swarm, numParticles, numUnknowns, particleLB, particleUB);
        }

        //Allocate Personal Particle Values
        double *particlePersonalBestValues = new double[numParticles];
        for (int i = 0; i < numParticles; i++)
        {
            particlePersonalBestValues[i] = numeric_limits<double>::max();
        }

        //This defaults to 0 which is what we want
        double *particlePersonalBestParticles = new double[numParticles * numUnknowns];
        double *costFunctionVals = new double[numParticles];
        double *particleVelocities = new double[numParticles * numUnknowns];
        particleBounds VelocityUB = particleUB - particleLB;
        particleBounds VelocityLB = particleLB - particleUB;
        double *velLBarray = VelocityLB.getBounds();
        double *velUBarray = VelocityUB.getBounds();
        double *posLBarray = particleLB.getBounds();
        double *posUBarray = particleUB.getBounds();

        //printParticleBests(particlePersonalBestValues, numParticles);
        //printSwarm(swarm, numParticles, numUnknowns);

        for (int iterationNum = 0; iterationNum < numIterations; iterationNum++)
        {
            /*
            ************Evaluate a single iteration within this loop************
            */

            if (displayIterationNum)
                cout << "***On iteration " << iterationNum + 1 << " out of " << numIterations << " ***\n";

            //We can parallelize here
            for (int particleNum = 0; particleNum < numParticles; particleNum++)
            {

                //Evaluate each particle here
                double deltaT1Particle = swarm[indexConversion(particleNum, 8, numUnknowns)];
                double xi0 = swarm[indexConversion(particleNum, 0, numUnknowns)];
                double xi1 = swarm[indexConversion(particleNum, 1, numUnknowns)];
                double xi2 = swarm[indexConversion(particleNum, 2, numUnknowns)];
                double xi3 = swarm[indexConversion(particleNum, 3, numUnknowns)];

                state_type inoutTarc1 = {tarc1.vrInitial, tarc1.vThetaInital, tarc1.rInitial, tarc1.xiInital};
                double t_start1 = 0.0, t_end1 = deltaT1Particle, dt = initialStepSize;
                //[ dense_output_detail_generation1

                /*
                @parameter absolute error tolerance
                @parameter relative error tolerance
            */

                dense_stepper_type tarc1Stepper = make_dense_output(absoluteTolerance, relTolerance, runge_kutta_dopri5<state_type>());

                if (outputParticleParams)
                {
                    cout << "xi0 is " << xi0 << endl;
                    cout << "xi1 is " << xi1 << endl;
                    cout << "xi2 is " << xi2 << endl;
                    cout << "xi3 is " << xi3 << endl;
                    cout << "Delta T1 particle is " << deltaT1Particle << endl
                         << endl;
                }

                thrustArc1EOM sys1 = thrustArc1EOM(ub, xi0, xi1, xi2, xi3, fuel);

                double vr1, vTheta1, r1, xi1PostTarc;
                try
                {
                    integrate_adaptive(tarc1Stepper, sys1, inoutTarc1, t_start1, t_end1, dt, timeStepObserver(dt));
                    vr1 = inoutTarc1[0];
                    vTheta1 = inoutTarc1[1];
                    r1 = inoutTarc1[2];
                    xi1PostTarc = inoutTarc1[3];
                }
                catch (MinStepSizeException &exception)
                {
                    vr1 = exception.getVr();
                    vTheta1 = exception.getvTheta();
                    r1 = exception.getR();
                    xi1PostTarc = exception.getXi();
                    if (outputWarnings)
                    {
                        cout << "Minimum step size reached in first integration" << endl;
                    }
                }

                double aCoast = ub * r1 / (2 * ub - r1 * (pow(vr1, 2) + pow(vTheta1, 2)));
                double eCoast = sqrt(1 - pow(r1, 2) * pow(vTheta1, 2) / (ub * aCoast));

                if (debug)
                {
                    cout << "vr1 is " << vr1 << endl;
                    cout << "vTheta1 is " << vTheta1 << endl;
                    cout << "r1 is " << r1 << endl;
                    cout << "xi1 is " << xi1PostTarc << endl;
                    cout << "aCoast is " << aCoast << endl;
                    cout << "e Coast is " << eCoast << endl
                         << endl;
                }

                if (aCoast > 0)
                {
                    double sinTrueAnamoly1 = vr1 / eCoast * sqrt(aCoast * (1 - pow(eCoast, 2)) / ub);
                    double cosTrueAnamoly1 = vTheta1 / eCoast * sqrt(aCoast * (1 - pow(eCoast, 2)) / ub) - 1 / eCoast;
                    double trueAnamoly1 = atan2(sinTrueAnamoly1, cosTrueAnamoly1);
                    if (trueAnamoly1 < 0)
                    {
                        trueAnamoly1 = trueAnamoly1 + 2 * M_PI;
                    }

                    double sinEccAnomaly1 = sin(trueAnamoly1) * sqrt(1 - pow(eCoast, 2)) / (1 + eCoast * cosTrueAnamoly1);
                    double cosEccAnomaly1 = (cos(trueAnamoly1) + eCoast) / (1 + eCoast * cosTrueAnamoly1);
                    double eccAnamoly1 = atan2(sinEccAnomaly1, cosEccAnomaly1);
                    if (eccAnamoly1 < 0)
                    {
                        eccAnamoly1 = eccAnamoly1 + 2 * M_PI;
                    }

                    double deltaE = swarm[indexConversion(particleNum, 9, numUnknowns)];

                    double eccAnomaly2 = eccAnamoly1 + deltaE;

                    double sinTrueAnamoly2 = sin(eccAnomaly2) * sqrt(1 - pow(eCoast, 2)) / (1 - eCoast * cos(eccAnomaly2));
                    double cosTrueAnamoly2 = (cos(eccAnomaly2) - eCoast) / (1 - eCoast * cos(eccAnomaly2));

                    double trueAnamoly2 = atan2(sinTrueAnamoly2, cosTrueAnamoly2);

                    if (trueAnamoly2 < 0)
                    {
                        trueAnamoly2 = trueAnamoly2 + 2 * M_PI;
                    }
                    double coastingTimeInterval = sqrt(pow(aCoast, 3) / ub) * (eccAnomaly2 - eccAnamoly1 - eCoast * (sin(eccAnomaly2) - sin(eccAnamoly1)));

                    //The first section is checked and working for good particles
                    if (debug)
                    {
                        cout << "True anamoly 1 is " << trueAnamoly1 << endl;
                        cout << "Eccentric anamoly 1 is " << eccAnamoly1 << endl;
                        cout << "Eccentric anamoly 2 is " << eccAnomaly2 << endl;
                        cout << "True anamoly 2 is " << trueAnamoly2 << endl;
                        cout << "Coasting time interval 2 is " << coastingTimeInterval << endl
                             << endl;
                    }

                    double vr2 = sqrt(ub / (aCoast * (1 - pow(eCoast, 2)))) * eCoast * sin(trueAnamoly2);
                    double vTheta2 = sqrt(ub / (aCoast * (1 - pow(eCoast, 2)))) * (1 + eCoast * cos(trueAnamoly2));
                    double r2 = aCoast * (1 - pow(eCoast, 2)) / (1 + eCoast * cos(trueAnamoly2));
                    double xi2PreTarc2 = xi1PostTarc + (trueAnamoly2 - trueAnamoly1);

                    //This is checked and working
                    if (debug)
                    {
                        cout << "vr2 is " << vr2 << endl;
                        cout << "vTheta2 is " << vTheta2 << endl;
                        cout << "r2 is " << r2 << endl;
                        cout << "xi2 is " << xi2PreTarc2 << endl
                             << endl;
                    }

                    double deltaT2Particle = swarm[indexConversion(particleNum, 10, numUnknowns)];
                    double v0 = swarm[indexConversion(particleNum, 4, numUnknowns)];
                    double v1 = swarm[indexConversion(particleNum, 5, numUnknowns)];
                    double v2 = swarm[indexConversion(particleNum, 6, numUnknowns)];
                    double v3 = swarm[indexConversion(particleNum, 7, numUnknowns)];

                    if (outputParticleParams)
                    {
                        cout << "v0 is " << v0 << endl;
                        cout << "v1 is " << v1 << endl;
                        cout << "v2 is " << v2 << endl;
                        cout << "v3 is " << v3 << endl;
                        cout << "deltaE is " << deltaE << endl;
                        cout << "deltaT2 is " << deltaT2Particle << endl
                             << endl;
                    }

                    thrustArcIc tarc2 = thrustArcIc(vr2, vTheta2, r2, xi2PreTarc2);
                    dense_stepper_type tarc2Stepper = make_dense_output(absoluteTolerance, relTolerance, runge_kutta_dopri5<state_type>());
                    double t_start2 = 0.0 + deltaT1Particle + coastingTimeInterval;
                    double t_end2 = deltaT2Particle + t_start2, dt = initialStepSize;
                    state_type inoutTarc2 = {tarc2.vrInitial, tarc2.vThetaInital, tarc2.rInitial, tarc2.xiInital};
                    thrustArc2EOM sys2 = thrustArc2EOM(ub, v0, v1, v2, v3, fuel, deltaT1Particle, coastingTimeInterval);

                    double vrFinal, vThetaFinal, rFinal, xiFinal;
                    try
                    {
                        integrate_adaptive(tarc2Stepper, sys2, inoutTarc2, t_start2, t_end2, dt, timeStepObserver(dt));
                        vrFinal = inoutTarc2[0];
                        vThetaFinal = inoutTarc2[1];
                        rFinal = inoutTarc2[2];
                        xiFinal = inoutTarc2[3];
                    }
                    catch (MinStepSizeException &exception)
                    {
                        vrFinal = exception.getVr();
                        vThetaFinal = exception.getvTheta();
                        rFinal = exception.getR();
                        xiFinal = exception.getXi();
                        if (outputWarnings)
                        {
                            cout << "Minimum step size reached in first integration" << endl;
                        }
                    }

                    if (debug)
                    {
                        cout << "vr Final is " << vrFinal << endl;
                        cout << "VTheta Final is " << vThetaFinal << endl;
                        cout << "rFinal is " << rFinal << endl;
                        cout << "xiFinal is " << xiFinal << endl
                             << endl;
                    }

                    double penalty = 0;

                    if (isnan(vrFinal) || isinf(vrFinal) || isnan(vThetaFinal) || isinf(vThetaFinal) || isnan(rFinal) || isinf(rFinal) || isnan(xiFinal) || isinf(xiFinal))
                    {
                        penalty += badResultPenalty;
                    }

                    if (abs(vrFinal) > penaltyValueCutoff)
                    {
                        penalty += abs(vrFinal) * penaltyCoefficient;
                    }

                    double d2 = vThetaFinal - sqrt(ub / R2);
                    if (abs(d2) > penaltyValueCutoff)
                    {
                        penalty += abs(d2) * penaltyCoefficient;
                    }

                    double d3 = rFinal - R2;
                    if (abs(d3) > penaltyValueCutoff)
                    {
                        penalty += abs(d3) * penaltyCoefficient;
                    }

                    if (debug)
                    {
                        cout << "vrError is " << vrFinal << endl;
                        cout << "vTheta Error is " << d2 << endl;
                        cout << "rFinal error is " << d3 << endl;
                        cout << "Penalty is " << penalty << endl
                             << endl;
                    }

                    double cost = deltaT1Particle + deltaT2Particle + penalty;
                    costFunctionVals[particleNum] = cost;
                }
                else
                {
                    costFunctionVals[particleNum] = badResultPenalty;
                }
                if (debug)
                {
                    cout << "Cost function value is " << costFunctionVals[particleNum] << endl
                         << endl;
                }
            }

            for (int particle = 0; particle < numParticles; particle++)
            {
                if (costFunctionVals[particle] < particlePersonalBestValues[particle])
                {
                    particlePersonalBestValues[particle] = costFunctionVals[particle];
                    setParticle(particlePersonalBestParticles, swarm, particle, numUnknowns);
                }
                if (costFunctionVals[particle] < globalBestValue)
                {
                    globalBestValue = costFunctionVals[particle];
                    setGlobalBestParticle(globalBestParticle, swarm, particle, numUnknowns);
                }
            }
            globalBestValuePerIteration[iterationNum] = globalBestValue;

            if (displayGlobalBestPerIteration)
            {
                cout << "The global best value for iteration number " << iterationNum + 1 << " is " << globalBestValue << endl;
            }

            //Now we update velocity

            double cI = (1 + randNumberPSO(0.0, 1.0)) / 2.0;
            double cC = 1.49445 * randNumberPSO(0.0, 1.0);
            double cS = 1.49445 * randNumberPSO(0.0, 1.0);

            if (displayPSOParams)
            {
                cout << "cI is " << cI << " : cC is " << cC << " : cS is " << cS << endl;
            }

            double cITerm, cCTerm, cSTerm, totalVelocityTerm;
            double *particleUBarray = particleUB.getBounds();
            double *particleLBarray = particleLB.getBounds();

            for (int i = 0; i < numParticles; i++)
            {
                for (int j = 0; j < numUnknowns; j++)
                {
                    int particleVariableIndex = indexConversion(i, j, numUnknowns);
                    cITerm = particleVelocities[particleVariableIndex] * cI;
                    cCTerm = cC * (particlePersonalBestParticles[particleVariableIndex] - swarm[particleVariableIndex]);
                    cSTerm = cS * (globalBestParticle[j] - swarm[particleVariableIndex]);
                    totalVelocityTerm = cITerm + cCTerm + cSTerm;

                    //Velocity bounds are correct
                    if (totalVelocityTerm < velLBarray[j])
                    {
                        totalVelocityTerm = velLBarray[j];
                    }
                    if (totalVelocityTerm > velUBarray[j])
                    {
                        totalVelocityTerm = velUBarray[j];
                    }
                    particleVelocities[particleVariableIndex] = totalVelocityTerm;
                }
            }

            //Update position here
            for (int i = 0; i < numParticles; i++)
            {
                for (int j = 0; j < numUnknowns; j++)
                {
                    //Position bounds are correct
                    int particleVariableIndex = indexConversion(i, j, numUnknowns);
                    swarm[particleVariableIndex] += particleVelocities[particleVariableIndex];
                    if (swarm[particleVariableIndex] < posLBarray[j])
                    {
                        swarm[particleVariableIndex] = posLBarray[j];
                    }
                    if (swarm[particleVariableIndex] > posUBarray[j])
                    {
                        swarm[particleVariableIndex] = posUBarray[j];
                    }
                }
            }
            if (displayCostFunctionAverages)
            {
                double avg = 0;
                for (int i = 0; i <= iterationNum; i++)
                {
                    avg += costFunctionVals[i];
                }
                avg /= static_cast<double>(iterationNum + 1);
                cout << "Cost function value average for iteration " << iterationNum + 1 << " is " << avg;
            }
        }

        cout << setprecision(16);
        printValues(globalBestParticle, numUnknowns, string("Global Best Particle is: "), numUnknowns);
        cout << "The global best value is: " << globalBestValue << endl;

        double wallFinal = getWallTime();
        double cpuFinal = getCPUTime();

        double wallElapsed = wallFinal - wall0;
        double cpuElapsed = cpuFinal - cpu0;

        cout << "Wall Time = " << wallElapsed << " seconds" << endl;
        cout << "CPU Time  = " << cpuElapsed << " seconds" << endl;

        writeToFile(resultsFile, wallElapsed, cpuElapsed, globalBestValue, globalBestParticle, numUnknowns, numParticles, numIterations);

        if (!evalTestParticle)
        {
            delete swarm;
        }
        delete globalBestParticle;
        delete particlePersonalBestValues;
        delete particlePersonalBestParticles;
        delete costFunctionVals;
        delete particleVelocities;
        delete globalBestValuePerIteration;
    }

    return 0;
}

void printSwarm(double *swarm, int numParticles, int numUnknowns)
{
    printValues(swarm, numParticles * numUnknowns, string("SWARM"), numUnknowns);
}

void printValues(double *values, int numToPrint, string description, int lineBreak = 6)
{
    string particleBestsString = "\n*************" + description + "************\n";
    int lineBreakCount = 0;
    cout << particleBestsString << endl;
    for (int i = 0; i < numToPrint; i++)
    {
        cout << values[i] << " ";
        lineBreakCount++;
        if (lineBreakCount % lineBreak == 0)
            cout << endl;
    }
    cout << particleBestsString << endl;
}

void allocateSwarm(double *swarm, int numParticles, int numUnknowns, particleBounds lowerBounds, particleBounds upperBounds)
{
    double *lb = lowerBounds.getBounds();
    double *ub = upperBounds.getBounds();
    for (int i = 0; i < numParticles; i++)
    {
        for (int j = 0; j < numUnknowns; j++)
        {
            //Rand should go between 0 and 1
            double rand = randNumberPSO(0.0, 1.0);
            double randBetweenBounds = (ub[j] - lb[j]) * rand + lb[j];
            swarm[indexConversion(i, j, numUnknowns)] = randBetweenBounds;
        }
    }
}

double randNumberPSO(double min, double max)
{
    double r3 = min + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX / (max - min)));
    return r3;
}

//(3,4) will return 4*11+3 = 44
int indexConversion(int x, int y, int arrWidth = numUnknowns)
{
    return x * arrWidth + y;
}

void setParticle(double *target, double *swarm, int numParticle, int numUnknowns)
{
    for (int i = 0; i < numUnknowns; i++)
    {
        int index = indexConversion(numParticle, i, numUnknowns);
        double value = swarm[index];
        target[index] = value;
    }
    return;
}

void setGlobalBestParticle(double *globalBestParticle, double *swarm, int numParticle, int numUnknowns)
{
    for (int i = 0; i < numUnknowns; i++)
    {
        globalBestParticle[i] = swarm[indexConversion(numParticle, i, numUnknowns)];
    }
    return;
}

double getWallTime()
{
    struct timeval time;
    if (gettimeofday(&time, NULL))
    {
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double getCPUTime()
{
    return (double)clock() / CLOCKS_PER_SEC;
}

void writeToFile(ofstream &fileStream, double wallElapsed, double cpuElapsed, double globalBestValue, double *globalBestParticle, int numUnknowns, int numParticles, int numIterations)
{
    fileStream << numParticles << "," << numIterations << ",";
    fileStream << wallElapsed << "," << cpuElapsed << "," << globalBestValue;
    for (int i = 0; i < numUnknowns; i++)
    {
        fileStream << "," << globalBestParticle[i];
    }
    fileStream << endl;
}

void writeHeadersToFile(ofstream &outFile){
    outFile << "Num Particles" << "," << "Num Iterations" << "," << "Wall elapsed (s)" << "," << "CPU Elapsed (s)" << ",";
    outFile << "JBest" << "," << "xi1" << "," << "xi2" << "," << "xi3" << "," << "xi4" << ",";
    outFile << "vi1" << "," << "vi2" << "," << "vi3" << "," << "vi4" << "," << "t1" << "," << "deltaE" << "," << "t2" << endl;
}



bool fileExists(string fileName)
{
    ifstream infile(fileName);
    return infile.good();
}