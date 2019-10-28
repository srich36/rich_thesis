vrInitial = 0; vrTerminal = 0;

R1 = 100; R2 = 200; ub = 500;
rInitial = R1; rFinal = R2;

vThetaInital = sqrt(ub/R1); vThetaFinal = sqrt(ub/R2);

spaceCraftMassInitial = 1000; thrustLevel = 500;
n0 = thrustLevel/spaceCraftMassInitial;


% Equation 57
aCoast = ub*r1/(2*ub-r1*(vr1^2+vTheta1^2));

%Equation 58
eCoast = sqrt(1-r1^2*vTheta1^2/(ub*aCoast));

%Equation 59
sinTrueAnamoly1 = vr1/eCoast*sqrt(aCoast*(1-eCoast)/ub);
cosTrueAnamoly1 = vTheta1/eCoast*sqrt(aCoast*(1-eCoast)/ub)-1/eCoast;

trueAnamoly1 = atan2(sinTrueAnamoly1, cosTrueAnamoly1);

% Equation 60
sinEccAnomaly1 = sinTrueAnamoly1*sqrt(1-eCoast^2)/(1+eCoast*cosTrueAnamoly1);
cosEccAnomaly1 = (cosTrueAnamoly1+eCoast)/(1+eCoast*cosTrueAnamoly1);

eccAnamoly1 = atan2(sinEccAnomaly1, cosEccAnomaly1);


deltaE = .5;

eccAnomaly2 = eccAnamoly1+deltaE;

%Equation 61
sinTrueAnamoly2 = sin(eccAnomaly2)*sqrt(1-eCoast^2)/(1-eCoast*cos(eccAnomaly2));
cosTrueAnamoly2 = (cos(eccAnomaly2)-eCoast)/(1-eCoast*cos(eccAnomaly2));

trueAnamoly2 = atan2(sinTrueAnamoly2, cosTrueAnamoly2);
%Equation 62
coastingTimeInterval = sqrt(aCoast^3/ub)*(eccAnomaly2-eccAnamoly1-eCoast*(sin(eccAnomaly2)-sinEccAnomaly1);


%{
    Inital Conditions for thrust arc 2
    (equations 63-66)
%}

%Equation 63
vr2 = sqrt(ub/(aCoast*(1-eCoast^2)))*eCoast*sinTrueAnamoly2;

%Equation 64
vtheta2 = sqrt(ub/(aCoast*(1-eCoast^2)))*(1+eCoast*cosTrueAnamoly2);

%Equation 65
r2 = aCoast*(1-eCoast^2)/(1+eCoast*cosTrueAnamoly2);

%Equation 66

xi2 = xi1 +(trueAnamoly2-trueAnamoly1);
