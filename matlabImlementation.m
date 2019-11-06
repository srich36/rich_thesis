clc; clear all; close;

penaltyCoefficient = 100;
vrInitial = 0; vrTerminal = 0;

xiInital = 0;
Beta = 2;
R1 = 1; R2 = Beta*R1; ub = 1;
rInitial = R1; rFinal = R2;

vThetaInital = sqrt(ub/R1); vThetaFinal = sqrt(ub/R2);

spaceCraftMassInitial = 1000; thrustLevel = 500;
n0 = thrustLevel/spaceCraftMassInitial;

icThurstArc1 = [ vrInitial; vThetaInital; rInitial; xiInital ];


%{
Parameter Bounds
%}

LBt1=0; LBdeltaE = 0; LBdeltat2=0; LBxi = -1; LBv=-1;
UBt1=3; UBdeltaE = 2*pi; UBdeltat2=3; UBxi = 1; UBv=1;

ParticleLB = [ LBxi, LBxi, LBxi, LBxi, LBv, LBv, LBv, LBv, LBt1, LBdeltaE, LBdeltat2 ];
ParticleUB = [ UBxi, UBxi, UBxi, UBxi, UBv, UBv, UBv, UBv, UBt1, UBdeltaE, UBdeltat2 ];

%{
    Assign Parameters
%}

numParticles = 110; numUnknowns = 11; numIterations = 100;

swarm = zeros(numParticles, numUnknowns);

for i=1:numUnknowns
    swarm(:, i) = ParticleLB(i)+rand(numParticles,1)*(ParticleUB(i)-ParticleLB(i));
end
pBest = zeros(numParticles, numUnknowns);
costFunctionVals = zeros(numParticles, 1);
particleCostFunctionBests = zeros(numParticles, 1);
velocities = zeros(numParticles, numUnknowns);

VelocityUB = ParticleUB-ParticleLB;
VelocityLB = -VelocityUB;

for i = 1:numParticles
    JBestVals(i) = inf;
end

globalBestValue = inf;
globalVestValuePerIteration = zeros(numIterations, 1);
globalBestParticle = zeros(1, numUnknowns);

for j=1:numIterations

    %Functions to evaluate each particle

    %Integrate the first round here
    for k=1:numParticles
        particle = swarm(k, :)
        deltaT1Particle = particle(9);
        tSpan = 0:.001:deltaT1Particle;
        tSpan = [0 deltaT1Particle]

        %Initial conditions are for vr, vtheta, r, psi
        iC1 = [ vrInitial; vThetaInital; rInitial; xiInital ];
        [t,thArc1Results] = ode45(@eomSolver1,tSpan, iC1);
        vr1 = thArc1Results(1);
        vTheta1 = thArc1Results(2);
        r1 = thArc1Results(3);
        xi1 = thArc1Results(4);


        % Equation 57
        aCoast = ub*r1/(2*ub-r1*(vr1^2+vTheta1^2));

        %Equation 58
        eCoast = sqrt(1-r1^2*vTheta1^2/(ub*aCoast));

        %Equation 59
        sinTrueAnamoly1 = vr1/eCoast*sqrt(aCoast*(1-eCoast)/ub);
        cosTrueAnamoly1 = vTheta1/eCoast*sqrt(aCoast*(1-eCoast)/ub)-1/eCoast;

        trueAnamoly1 = atan2(sinTrueAnamoly1, cosTrueAnamoly1);
        if trueAnamoly1 < 0
        trueAnamoly1+=pi;
        end

        % Equation 60
        sinEccAnomaly1 = sinTrueAnamoly1*sqrt(1-eCoast^2)/(1+eCoast*cosTrueAnamoly1);
        cosEccAnomaly1 = (cosTrueAnamoly1+eCoast)/(1+eCoast*cosTrueAnamoly1);

        eccAnamoly1 = atan2(sinEccAnomaly1, cosEccAnomaly1);
        if eccAnamoly1 < 0
        eccAnamoly1+=pi;
        end

        deltaE = .5;

        eccAnomaly2 = eccAnamoly1+deltaE;

        %Equation 61
        sinTrueAnamoly2 = sin(eccAnomaly2)*sqrt(1-eCoast^2)/(1-eCoast*cos(eccAnomaly2));
        cosTrueAnamoly2 = (cos(eccAnomaly2)-eCoast)/(1-eCoast*cos(eccAnomaly2));

        trueAnamoly2 = atan2(sinTrueAnamoly2, cosTrueAnamoly2);
        if trueAnamoly2 < 0
        trueAnamoly2+=pi;
        end
        %Equation 62
        coastingTimeInterval = sqrt(aCoast^3/ub)*(eccAnomaly2-eccAnamoly1-eCoast*(sin(eccAnomaly2)-sinEccAnomaly1));

        %{
        %    Inital Conditions for thrust arc 2
        %    (equations 63-66)
        %}

        %Equation 63
        vr2 = sqrt(ub/(aCoast*(1-eCoast^2)))*eCoast*sinTrueAnamoly2;

        %Equation 64
        vtheta2 = sqrt(ub/(aCoast*(1-eCoast^2)))*(1+eCoast*cosTrueAnamoly2);

        %Equation 65
        r2 = aCoast*(1-eCoast^2)/(1+eCoast*cosTrueAnamoly2);

        %Equation 66

        xi2 = xi1 + (trueAnamoly2-trueAnamoly1);

        deltaT2Particle = particle(11);
        iC2 = [ vr2; vtheta2; r2; xi2 ];
        [t,thArc2Results] = ode45(@eomSolver1,tSpan, iC1);
        vrFinal = thArc2Results(1);
        vThetaFinal = thArc2Results(2);
        rFinal = thArc2Results(3);
        xiFinal = thArc2Results(4);

        penalty = 0;
        if vrFinal > 10e-3
            penalty+=vrFinal*penaltyCoefficient;
        end

        d2 = vThetaFinal-sqrt(ub/R2);
        if d2 > 10e-3
            penalty+=d2*penaltyCoefficient;
        end

        d3 = rFinal-R2;
        if d3 > 10e-3
            penalty+=d3*penaltyCoefficient;
        end

        costFunctionVals(k) = deltaT1Particle+deltaT2Particle+penalty;

    end


    %costFunctionVals = evaluateUnknowns(swarm, costFunctionVals, numParticles);

    % Setting local and global best particles and values
    for i = 1:numParticles

        %Checking for local best particle and value
        if costFunctionVals(i) < particleCostFunctionBests(i)
            pBest(i, :) = swarm(i, :);
            particleCostFunctionBests(i) = costFunctionVals(i);
        end

        %Checking for global best particle and value
        if costFunctionVals(i) < GlobalBest
            globalBestValue = costFunctionVals(i);
            globalBestParticle = swarm(i, :)
        end

        %Update particles
        cI = (1+rand)/2;
        cC = 1.49445*rand;
        cS = 1.49445*rand;

        % Update particle velocity
        for i = 1:numParticles

            cITerm = cI*velocities(i, :);
            cCTerm = cC*(particleCostFunctionBests(i, :) - swarm(i, :));
            cSTerm = cS*(globalBestParticle-P(i, :));

            % Need to check to make sure that this bounded value works with matrices
            velocities(i, :) = boundValue(cITerm+cCTerm+cSTerm, VelocityLB, VelocityUB);


            swarm(i, :) = swarm(i, :) + velocities(i, :);
            for k=1:numUnknowns
                swarm(i, k) = boundValue(swarm(i, k), ParticleLB(k), ParticleUB(k));
            end

        end

    end

    globalVestValuePerIteration(j) = globalBestValue;
end

save 'PSO_results' globalVestValuePerIteration, numIterations



