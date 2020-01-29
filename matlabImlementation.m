clc; clear all; close;
%warning("on");
warning("off");

%%TOLERANCES AT 1E-9 AND 1E-9%%
absTol =  1e-9;
relTol = 1e-9;
odeOpts = odeset('RelTol',relTol,'AbsTol', absTol);

testParticle = [ 0.1793388403978833	-0.8666167293262341	0.9191858318050022	-0.03751637123140417	0.01516904721637628	0.3782818408262493	-0.467870454436066	-1	0.6724822850591763	2.821471317903069	0.4106645490059094 ];
%Good one
%testParticle = [ -0.01688571838897738	0.09587006665201843	-0.7129715153615238	0.8309410641294435	0.1915107672508798	-0.7028940222556719	-0.4393799017564173	0.6916400209722409	0.6702618694967299	2.825282566477961	0.4135842958459761 ] ;

evalTestParticle = true;

cyclesBeforeHydrateAgain = 10;
cyclesBeforeHydrateAgainValue = 10;
resetVelocityTerm = [ 0 0 0 0 0 0 0 0 0 0 0 ];

debug = false;
debugPSO = false;
rehydrate = false;
graphValuesOverIterations = false;
%colorValues = [ 'bo-', 'ro-', 'yo-', 'go-', 'co-' ];

averageJThresholdPercent = .1;
numIterBeforeHydration = 50;
numPrevIterationsAvgJ = 10;

rehydratePercentage = .33;
displayCostFunctionValue = false;
displayAverageCostValuePerIteration = false;
displayGlobalBestPerIteration = true;
displayIterationNum = false;
displayExecutionNum = true;
displayPsoIterationNumber = true;
displayGlobalBestPerPSOIteration = true;
plotThrustArcs = true;

badResultPenalty = 100000;

penaltyCoefficient = 100;
penaltyValueCutoff = .001;

vrInitial = 0; vrTerminal = 0;
xiInital = 0;
Beta = 2;
R1 = 1; R2 = Beta*R1; ub = 1;
rInitial = R1; rFinal = R2;

vThetaInital = sqrt(ub/R1); vThetaFinal = sqrt(ub/R2);

c = .5;
n0 = .2;

icThurstArc1 = [ vrInitial; vThetaInital; rInitial; xiInital ];


%{
Parameter Bounds
%}

LBt1=0.001; LBdeltaE = 0; LBdeltat2=0.01; LBxi = -1; LBv=-1;
UBt1=3; UBdeltaE = 2*pi; UBdeltat2=3; UBxi = 1; UBv=1;

ParticleLB = [ LBxi, LBxi, LBxi, LBxi, LBv, LBv, LBv, LBv, LBt1, LBdeltaE, LBdeltat2 ];
ParticleUB = [ UBxi, UBxi, UBxi, UBxi, UBv, UBv, UBv, UBv, UBt1, UBdeltaE, UBdeltat2 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Particle Params     %%
numExecutions=1; numPsoIterations=1;
numParticles = 110; numUnknowns = 11; numIterations = 1000;
%%    Particle Params     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

thurstArcFigFilename = 'figures/thrustArcs.fig';

if plotThrustArcs
    if isfile(thurstArcFigFilename)
        openfig('figures/thrustArcs.fig'); 
    end
end

rehydrateString = 'rehydrate';

mkdir './results'


betaDirString = strcat('beta',num2str(Beta));
mkdir(fullfile('.', 'results', betaDirString));
if rehydrate
    mkdir(fullfile('.', 'results', betaDirString, rehydrateString));
end
fileNameString = strcat('LoopNum',num2str(numPsoIterations), 'pNum', num2str(numParticles), 'Inum', num2str(numIterations),'r', num2str(rehydrate), '.csv');


if ~rehydrate
    resultsFileName=fullfile('.','results', betaDirString, fileNameString);
else 
    resultsFileName=fullfile('.','results', betaDirString, rehydrateString, fileNameString);
end

if ~isfile(resultsFileName)
    cHeader = {'numParticles', 'numIterations', 'Jbest' 'Time Elapsed' 'xi1' 'xi2' 'xi3' 'xi4' 'vi1' 'vi2' 'vi3' 'vi4' 't1' 'deltaE' 't2' }; %dummy header
    commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
    commaHeader = commaHeader(:)';
    textHeader = cell2mat(commaHeader); %cHeader in text with commas

    fid = fopen(resultsFileName,'a');
    fprintf(fid,'%s\n',textHeader);
    fclose(fid);
end


for executionNum=1:numExecutions
    globalBestValuePerExecution = inf;
    globalBestParticlePerExecution = zeros(1, numUnknowns);

    
    if displayExecutionNum
        fprintf("***On execution %f out of %f***\n", executionNum, numExecutions);
    end
    %Start timer here so it times each execution
    tic
    
    for psoIteration=1:numPsoIterations 
        %{
            Assign Parameters
        %}
        
        if displayPsoIterationNumber
            fprintf("***On pso Iteration number %f out of %f***\n", psoIteration, numPsoIterations);            
        end
        swarm = zeros(numParticles, numUnknowns);
        cCTerm = zeros(numParticles, numUnknowns);
        cITerm = zeros(numParticles, numUnknowns);
        cSTerm = zeros(numParticles, numUnknowns);


        for i=1:numUnknowns
            swarm(:, i) = ParticleLB(i)+rand(numParticles,1)*(ParticleUB(i)-ParticleLB(i));
        end
        
        if psoIteration > 1
            for i=1:numParticles/5
                swarm(i,:) = globalBestParticlePerExecution;
                
            end
        end
        particlePersonalBest = zeros(numParticles, numUnknowns);
        costFunctionVals = zeros(numParticles, 1);
        particleCostFunctionBests = zeros(numParticles, 1);
        velocities = zeros(numParticles, numUnknowns);

        VelocityUB = ParticleUB-ParticleLB;
        VelocityLB = -VelocityUB;

        if evalTestParticle
            swarm = testParticle;
            numParticles = 1;
            numIterations = 1;
        end

        for i = 1:numParticles
            particleCostFunctionBests(i) = inf;
        end

        globalBestValue = inf;
        globalBestValuePerIteration = zeros(numIterations, 1);
        globalBestParticle = zeros(1, numUnknowns);

        for j=1:numIterations
            costValueSum = 0;
            
            if displayIterationNum
                fprintf("***On iteration %f out of %f***\n",j, numIterations);
            end
            
            Jsum = 0;
            if j > numIterBeforeHydration
                if j > numPrevIterationsAvgJ
                    prevJToCompare = globalBestValuePerIteration(j-numPrevIterationsAvgJ);
                    for jsumIter=flip(1:numPrevIterationsAvgJ)
                        Jsum = Jsum + globalBestValuePerIteration(j-jsumIter);
                        
                    end
                    %Since you are comparing the two previous values, the
                    %number of averages you are summing is the n previous
                    %values-1
                    avgJ = Jsum/numPrevIterationsAvgJ;
                    averageJChangePercent = 100*(prevJToCompare-avgJ)/prevJToCompare;
                else
                    averageJChangePercent = 100;
                end
            else
                averageJChangePercent = 100;
                
            end
            if rehydrate && averageJChangePercent < averageJThresholdPercent && ~cyclesBeforeHydrateAgain
                for rehydratedParticle = 1:floor(numParticles*rehydratePercentage)
                    for i=1:numUnknowns
                        swarm(rehydratedParticle, i) = ParticleLB(i)+rand(1,1)*(ParticleUB(i)-ParticleLB(i));
                        velocities(rehydratedParticle,:) = resetVelocityTerm;
                    end
                end
                disp("Rehydrating");
                cyclesBeforeHydrateAgain = cyclesBeforeHydrateAgainValue;
            end
            
            
            %Functions to evaluate each particle
            %Integrate the first round here
            for k=1:numParticles
                particle = swarm(k, :);
                deltaT1Particle = particle(9);
                %tSpan = 0:.001:deltaT1Particle;
                tSpan = [0 deltaT1Particle];
                xi0 = particle(1);
                %Initial conditions are for vr, vtheta, r, psi
                [t,y] = ode45(@(t,y) eomSolver1(t,y,ub, c, n0, xi0, particle(2), particle(3), particle(4)),tSpan, icThurstArc1, odeOpts);
                vr1 = y(end, 1);
                vTheta1 = y(end, 2);
                r1 = y(end, 3);
                xi1 = y(end, 4);

                % Equation 57 - checked
                aCoast = ub*r1/(2*ub-r1*(vr1^2+vTheta1^2));
                %Equation 58 - checked
                eCoast = sqrt(1-r1^2*vTheta1^2/(ub*aCoast));

                if debug
                    fprintf("Vr1 is: %f\n", vr1);
                    fprintf("VTheta1 is: %f\n", vTheta1);
                    fprintf("r1 is: %f\n", r1);
                    fprintf("xi1 is: %f\n\n", xi1);
                    fprintf("a coast is: %f\n", aCoast);
                    fprintf("e Coast is %f\n", eCoast);
                end

                if aCoast > 0
                    %Equation 59
                    sinTrueAnamoly1 = vr1/eCoast*sqrt(aCoast*(1-eCoast^2)/ub); % checked 
                    cosTrueAnamoly1 = vTheta1/eCoast*sqrt(aCoast*(1-eCoast^2)/ub)-1/eCoast; %checked

                    trueAnamoly1 = atan2(sinTrueAnamoly1, cosTrueAnamoly1);
                    if trueAnamoly1 < 0
                    trueAnamoly1=trueAnamoly1+2*pi;
                    end


                    % Equation 60
                    sinEccAnomaly1 = sin(trueAnamoly1)*sqrt(1-eCoast^2)/(1+eCoast*cosTrueAnamoly1);
                    cosEccAnomaly1 = (cos(trueAnamoly1)+eCoast)/(1+eCoast*cosTrueAnamoly1);
                    eccAnamoly1 = atan2(sinEccAnomaly1, cosEccAnomaly1);
                    if eccAnamoly1 < 0
                    eccAnamoly1=eccAnamoly1+2*pi;
                    end


                    deltaE = particle(10);
                    eccAnomaly2 = eccAnamoly1+deltaE;

                    %Equation 61
                    sinTrueAnamoly2 = sin(eccAnomaly2)*sqrt(1-eCoast^2)/(1-eCoast*cos(eccAnomaly2));
                    cosTrueAnamoly2 = (cos(eccAnomaly2)-eCoast)/(1-eCoast*cos(eccAnomaly2));

                    trueAnamoly2 = atan2(sinTrueAnamoly2, cosTrueAnamoly2);

                    if trueAnamoly2 < 0
                    trueAnamoly2=trueAnamoly2+2*pi;
                    end

                    %Equation 62
                    coastingTimeInterval = sqrt(aCoast^3/ub)*(eccAnomaly2-eccAnamoly1-eCoast*(sin(eccAnomaly2)-sin(eccAnamoly1)));

                    if debug
                        fprintf("True anamoly 1 is %f\n", trueAnamoly1);
                        fprintf("Eccentric anamoly 1 is %f\n\n", eccAnamoly1);
                        fprintf("Eccentric anamoly 2 is %f\n", eccAnomaly2);
                        fprintf("True anamoly 2 is %f\n", trueAnamoly2);
                        fprintf("Coasting time interval is %f\n\n", coastingTimeInterval);
                    end

                    %{
                    %    Inital Conditions for thrust arc 2
                    %    (equations 63-66)
                    %}

                    %Equation 63
                    vr2 = sqrt(ub/(aCoast*(1-eCoast^2)))*eCoast*sin(trueAnamoly2);

                    %Equation 64
                    vtheta2 = sqrt(ub/(aCoast*(1-eCoast^2)))*(1+eCoast*cos(trueAnamoly2));

                    %Equation 65
                    r2 = aCoast*(1-eCoast^2)/(1+eCoast*cos(trueAnamoly2));

                    %Equation 66

                    xi2 = xi1 + (trueAnamoly2-trueAnamoly1);

                    if debug
                        fprintf("vr2 is %f\n", vr2);
                        fprintf("vtheta2 is %f\n", vtheta2);
                        fprintf("r2 is %f\n", r2);
                        fprintf("xi2 is %f\n\n", xi2);
                    end

                    deltaT2Particle = particle(11);
                    iC2 = [ vr2; vtheta2; r2; xi2 ];
                    t2 = deltaT1Particle+coastingTimeInterval;
                    tSpan2 = [ t2 t2+deltaT2Particle ];
                    [t,thArc2Results] = ode45(@(t,y) eomSolver2(t,y,ub, c, n0, deltaT1Particle, t2, particle(5), particle(6), particle(7), particle(8)),tSpan2, iC2, odeOpts);

                    vrFinal = thArc2Results(end, 1);
                    vThetaFinal = thArc2Results(end, 2);
                    rFinal = thArc2Results(end, 3);
                    xiFinal = thArc2Results(end, 4);
                    if debug
                        fprintf("vrFinal is %f\n", vrFinal);
                        fprintf("vThetaFinal is %f\n", vThetaFinal);
                        fprintf("rFinal is %f\n", rFinal);
                        fprintf("xiFinal is %f\n\n", xiFinal);
                    end

                    penalty = 0;

                    if isnan(vrFinal) || isnan(vThetaFinal) || isnan(rFinal) || isnan(xiFinal)
                       penalty = penalty+badResultPenalty;
                    end

                    if isinf(abs(vrFinal)) || isinf(abs(vThetaFinal)) || isinf(abs(rFinal)) || isinf(abs(xiFinal))
                       penalty = penalty+badResultPenalty;
                    end

                    if abs(vrFinal) > penaltyValueCutoff
                        penalty=penalty+abs(vrFinal)*penaltyCoefficient;
                    end

                    d2 = vThetaFinal-sqrt(ub/R2);
                    if abs(d2) > penaltyValueCutoff
                        penalty=penalty+abs(d2)*penaltyCoefficient;
                    end

                    d3 = rFinal-R2;
                    if abs(d3) > penaltyValueCutoff
                        penalty=penalty+abs(d3)*penaltyCoefficient;
                    end
                    if evalTestParticle
                        fprintf("vrFinal is %f\n", vrFinal)
                        fprintf("vTheta error value is %f\n", d2);
                        fprintf("rFinal error value is %f\n", d3);
                    end
                    cost = deltaT1Particle+deltaT2Particle+penalty;
                    costFunctionVals(k) = cost;
                else
                    costFunctionVals(k) = badResultPenalty;
                end
                if displayCostFunctionValue
                    fprintf("Cost function value is %f\n\n", costFunctionVals(k));
                end
                costValueSum=costValueSum+costFunctionVals(k);
            end

            costValueAverage = costValueSum/numParticles;
            if displayAverageCostValuePerIteration
                fprintf("For iteration %f out of %f the average cost function value is %f\n\n", j, numIterations, costValueAverage);
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
                    if globalBestValue < globalBestValuePerExecution
                       globalBestValuePerExecution = globalBestValue;
                       globalBestParticlePerExecution = globalBestParticle;
                    end
                end

            end

            % Update particle velocity
            for i = 1:numParticles
                %Update particles
                cI = (1+rand)/2;
                cC = 1.49445*rand;
                cS = 1.49445*rand;

                %Particle Velocity term
                cITerm(i,:) = cI*velocities(i, :);

                %Deviation from Personal Best term
                cCTerm(i,:) = cC*(particlePersonalBest(i, :) - swarm(i, :));

                %Deviation from Global Best Term
                cSTerm(i,:) = cS*(globalBestParticle-swarm(i, :));

                totalVelocityTerm = cSTerm(i,:)+cCTerm(i,:)+cITerm(i,:);
                velocities(i,:) = totalVelocityTerm;
                for k=1:numUnknowns
                    if velocities(i,k) > VelocityUB(k)
                        velocities(i,k) = VelocityUB(k);
                    end

                    if velocities(i,k) < VelocityLB(k)
                        velocities(i,k) = VelocityLB(k);
                    end
                end

                swarm(i, :) = swarm(i, :) + velocities(i, :);
                for k=1:numUnknowns
                    if swarm(i, k) < ParticleLB(k)
                       swarm(i,k) = ParticleLB(k); 
                    end

                    if swarm(i, k) > ParticleUB(k)
                       swarm(i,k) = ParticleUB(k); 
                    end
                end

            end

            if debugPSO
                fprintf("cI term is %f\n", cITerm);
                fprintf("cC term is %f\n", cCTerm);
                fprintf("cS term is %f\n", cSTerm);
            end

            globalBestValuePerIteration(j) = globalBestValue;
            if displayGlobalBestPerIteration
                fprintf("Global best for iteration %f is %f\n", j, globalBestValue);
            end
            
            if cyclesBeforeHydrateAgain > 0
                cyclesBeforeHydrateAgain = cyclesBeforeHydrateAgain - 1;
            end
        end
        
        if displayGlobalBestPerPSOIteration
            fprintf("Global best for PSO iteration %f is %f\n", psoIteration, globalBestValue);
        end
        
        if plotThrustArcs
            
            p = aCoast*(1-eCoast^2);
            xiCoast = xi1:.01:xiFinal;
            trueAnamolyCoast = trueAnamoly1:.01:trueAnamoly2;


            rCoastIndex = 1;
            for trueAnamoly = trueAnamolyCoast

                rCoast(rCoastIndex) = aCoast*(1-eCoast^2)/(1+eCoast*cos(trueAnamoly));
                rCoastIndex=rCoastIndex+1;
            end
            
            xiCoast = linspace(xi1, xiFinal, rCoastIndex-1);
            
            xiValues = 0:.01:2*pi;
            rIndex = 1;
            
            
            for xi=xiValues
                rInitial(rIndex) = 1;
                rFinal(rIndex) = Beta;
                rIndex = rIndex+1;
            end
            
            rFirstArc = y(:,3)';
            rLastArc = thArc2Results(:,3)';
            
            thetaFirstArc = y(:,4)';
            thetaLastArc = thArc2Results(:,4)';
            rValues = horzcat(rFirstArc, rCoast, rLastArc);
            thetaValues = horzcat(thetaFirstArc, xiCoast, thetaLastArc);
            
            polarplot(xiValues, rInitial);
            hold on;
            polarplot(thetaValues, rValues);
            hold on;
            polarplot(xiValues, rFinal);
            title("Optimal transfer trajectory (Beta=2)");
            legend("Initial orbit", "Transfer Trajectory", "Final orbit");
            
        end
    end
    executionTime = toc;
    
    if graphValuesOverIterations
        plot(1:numIterations, globalBestValuePerIteration, '-o');
        title("JBest Value Per Iteration - 110 Particles");
        xlabel("Iteration number");
        ylabel("Global JBest Value");
        hold on;
        legend('Run 1', 'Run 2', 'Run 3', 'Run 4', 'Run 5');
    end
    disp(globalBestParticlePerExecution);
    csvData = horzcat(numParticles, numIterations, globalBestValuePerExecution, executionTime, globalBestParticlePerExecution);
    dlmwrite(resultsFileName,csvData,'-append', 'precision', 16);

    
end

save 'PSO_results_matlab'



fprintf("\nDelta t1 for best particle is %f\n", globalBestParticlePerExecution(9));
fprintf("Delta E for best particle is %f\n", globalBestParticlePerExecution(10));
fprintf("Delta t2 for best particle is %f\n\n", globalBestParticlePerExecution(11));
disp(globalBestParticlePerExecution);



