clc; clear all; close;
warning("off");

debug = false;
displayCostFunctionValue = false;
displayAverageCostValuePerIteration = true;

penaltyCoefficient = 100;
vrInitial = 0; vrTerminal = 0;

xiInital = 0;
Beta = 2;
R1 = 1; R2 = Beta*R1; ub = 1;
rInitial = R1; rFinal = R2;

vThetaInital = sqrt(ub/R1); vThetaFinal = sqrt(ub/R2);

%spaceCraftMassInitial = 1000; thrustLevel = 500;
c = .5;
n0 = .2;

icThurstArc1 = [ vrInitial; vThetaInital; rInitial; xiInital ];


%{
Parameter Bounds
%}

LBt1=0.001; LBdeltaE = 0; LBdeltat2=0.001; LBxi = -1; LBv=-1;
UBt1=3; UBdeltaE = 2*pi; UBdeltat2=3; UBxi = 1; UBv=1;

ParticleLB = [ LBxi, LBxi, LBxi, LBxi, LBv, LBv, LBv, LBv, LBt1, LBdeltaE, LBdeltat2 ];
ParticleUB = [ UBxi, UBxi, UBxi, UBxi, UBv, UBv, UBv, UBv, UBt1, UBdeltaE, UBdeltat2 ];

%{
    Assign Parameters
%}

numParticles = 30; numUnknowns = 11; numIterations = 10;

swarm = zeros(numParticles, numUnknowns);

for i=1:numUnknowns
    swarm(:, i) = ParticleLB(i)+rand(numParticles,1)*(ParticleUB(i)-ParticleLB(i));
end
particlePersonalBest = zeros(numParticles, numUnknowns);
costFunctionVals = zeros(numParticles, 1);
particleCostFunctionBests = zeros(numParticles, 1);
velocities = zeros(numParticles, numUnknowns);

VelocityUB = ParticleUB-ParticleLB;
VelocityLB = -VelocityUB;

for i = 1:numParticles
    particleCostFunctionBests(i) = inf;
end

globalBestValue = inf;
globalBestValuePerIteration = zeros(numIterations, 1);
globalBestParticle = zeros(1, numUnknowns);

for j=1:numIterations
    costValueSum = 0;
    printf("***On iteration %f out of %f***\n",j, numIterations)
    %Functions to evaluate each particle

    %Integrate the first round here
    for k=1:numParticles
        particle = swarm(k, :);
        deltaT1Particle = particle(9);
        %tSpan = 0:.001:deltaT1Particle;
        tSpan = [0 deltaT1Particle];
        xi0 = particle(1);
        %Initial conditions are for vr, vtheta, r, psi
        [t,y] = ode45(@(t,y) eomSolver1(t,y,ub, c, n0, xi0, particle(2), particle(3), particle(4)),tSpan, icThurstArc1);
        vr1 = y(end, 1);
        vTheta1 = y(end, 2);
        r1 = y(end, 3);
        xi1 = y(end, 4);

        % Equation 57
        aCoast = ub*r1/(2*ub-r1*(vr1^2+vTheta1^2));
        %Equation 58
        eCoast = sqrt(1-r1^2*vTheta1^2/(ub*aCoast));
        if debug
            printf("Vr1 is: %f\n", vr1);
            printf("VTheta1 is: %f\n", vTheta1);
            printf("r1 is: %f\n", r1);
            printf("xi1 is: %f\n\n", xi1);
            printf("a coast is: %f\n", aCoast);
            printf("e Coast is %f\n", eCoast);
        end

        if aCoast > 0
            %Equation 59
            sinTrueAnamoly1 = vr1/eCoast*sqrt(aCoast*(1-eCoast)/ub);
            cosTrueAnamoly1 = vTheta1/eCoast*sqrt(aCoast*(1-eCoast)/ub)-1/eCoast;

            trueAnamoly1 = atan2(sinTrueAnamoly1, cosTrueAnamoly1);
            if trueAnamoly1 < 0
            trueAnamoly1+=2*pi;
            end

            % Equation 60
            sinEccAnomaly1 = sinTrueAnamoly1*sqrt(1-eCoast^2)/(1+eCoast*cosTrueAnamoly1);
            cosEccAnomaly1 = (cosTrueAnamoly1+eCoast)/(1+eCoast*cosTrueAnamoly1);
            eccAnamoly1 = atan2(sinEccAnomaly1, cosEccAnomaly1);
            if eccAnamoly1 < 0
            eccAnamoly1+=2*pi;
            end


            deltaE = particle(10);

            eccAnomaly2 = eccAnamoly1+deltaE;

            %Equation 61
            sinTrueAnamoly2 = sin(eccAnomaly2)*sqrt(1-eCoast^2)/(1-eCoast*cos(eccAnomaly2));
            cosTrueAnamoly2 = (cos(eccAnomaly2)-eCoast)/(1-eCoast*cos(eccAnomaly2));

            trueAnamoly2 = atan2(sinTrueAnamoly2, cosTrueAnamoly2);

            if trueAnamoly2 < 0
            trueAnamoly2+=2*pi;
            end

            %Equation 62
            coastingTimeInterval = sqrt(aCoast^3/ub)*(eccAnomaly2-eccAnamoly1-eCoast*(sin(eccAnomaly2)-sinEccAnomaly1));

            if debug
                printf("True anamoly 1 is %f\n", trueAnamoly1);
                printf("Eccentric anamoly 1 is %f\n\n", eccAnamoly1);
                printf("Eccentric anamoly 2 is %f\n", eccAnomaly2);
                printf("True anamoly 2 is %f\n", trueAnamoly2);
                printf("Coasting time interval is %f\n\n", coastingTimeInterval);
            end

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

            if debug
                printf("vr2 is %f\n", vr2);
                printf("vtheta2 is %f\n", vtheta2);
                printf("r2 is %f\n", r2);
                printf("xi2 is %f\n\n", xi2);
            end

            deltaT2Particle = particle(11);
            iC2 = [ vr2; vtheta2; r2; xi2 ];
            t2 = deltaT1Particle+coastingTimeInterval;
            [t,thArc2Results] = ode45(@(t,y) eomSolver2(t,y,ub, c, n0, deltaT1Particle, t2, particle(5), particle(6), particle(7), particle(8)),tSpan, iC2);
            vrFinal = thArc2Results(end, 1);
            vThetaFinal = thArc2Results(end, 2);
            rFinal = thArc2Results(end, 3);
            xiFinal = thArc2Results(end, 4);
            if debug
                printf("vrFinal is %f\n", vrFinal);
                printf("vThetaFinal is %f\n", vThetaFinal);
                printf("rFinal is %f\n", rFinal);
                printf("xiFinal is %f\n\n", xiFinal);
            end

            penalty = 0;
            if abs(vrFinal) > 10e-3
                penalty+=abs(vrFinal)*penaltyCoefficient;
            end

            d2 = vThetaFinal-sqrt(ub/R2);
            if abs(d2) > 10e-3
                penalty+=abs(d2)*penaltyCoefficient;
            end

            d3 = rFinal-R2;
            if abs(d3) > 10e-3
                penalty+=abs(d3)*penaltyCoefficient;
            end
            cost = deltaT1Particle+deltaT2Particle+penalty;
            costFunctionVals(k) = cost;
        else
            costFunctionVals(k) = 100000;
        end
        if displayCostFunctionValue
            printf("Cost function value is %f\n\n", costFunctionVals(k));
        end
        costValueSum+=costFunctionVals(k);
    end

    costValueAverage = costValueSum/numParticles;
    if displayAverageCostValuePerIteration
        printf("For iteration %f out of %f the average cost function value is %f\n\n", j, numIterations, costValueAverage);
    end
    % Setting local and global best particles and values
    for numPart = 1:numParticles

        %Checking for local best particle and value
        if costFunctionVals(numPart) < particleCostFunctionBests(numPart)
            particlePersonalBest(numPart, :) = swarm(numPart, :);
            particleCostFunctionBests(numPart) = costFunctionVals(numPart);
        end

        %Checking for global best particle and value
        if costFunctionVals(numPart) < globalBestValue
            globalBestValue = costFunctionVals(numPart);
            globalBestParticle = swarm(numPart, :);
        end

    end

    % Update particle velocity
    for i = 1:numParticles
        %Update particles
        cI = (1+rand)/2;
        cC = 1.49445*rand;
        cS = 1.49445*rand;
        cITerm = cI*velocities(i, :);
        cCTerm = cC*(particlePersonalBest(i, :) - swarm(i, :));
        disp(cCTerm);
        cSTerm = cS*(globalBestParticle-swarm(i, :));

        if debug
            printf("cI term is %f\n", cITerm);
            printf("cC term is %f\n", cCTerm);
            printf("cS term is %f\n", cSTerm);
        end
        % Need to check to make sure that this bounded value works with matrices
        velocities(i, :) = boundValue(cITerm+cCTerm+cSTerm, VelocityLB, VelocityUB);
        swarm(i, :) = swarm(i, :) + velocities(i, :);
        for k=1:numUnknowns
            swarm(i, k) = boundValue(swarm(i, k), ParticleLB(k), ParticleUB(k));
        end

    end

    globalBestValuePerIteration(j) = globalBestValue;
    printf("Global best for iteration %f is %f\n", j, globalBestValue);
end

save 'PSO_results'
disp(globalBestParticle);


