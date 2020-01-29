particlesB2 = [ -0.1901198028497241	1	-0.04834424815088992	-1	-0.0003163896310721081	-0.3263155003733011	0.8848431828601635	-0.09326100639757344	0.6709061443951977	2.76405176917629	0.4132477279882463;
    0.16869	-0.7152	0.13536	1	0.13758	-0.38264	-1	0.28981	0.67163	2.7975	0.41235;
    0.024413	0.24366	0.032103	-0.43868	0.1351	-0.74577	0.41254	-0.3803	0.67002	2.7597	0.41335 ];



particlesB4 = [-0.083403 0.66131 -0.44825	-0.081968	0.12607	-0.34315	-0.13379	-0.57554	1.046	2.8183	0.43996;
   -0.14636	0.10731	0.28391	-0.067937	-0.039543	0.50864	0.32392	0.19839	1.0434	2.9085	0.44457;
    -0.23143 0.58274	0.18055	-0.42708	0.24862	-0.9965	0.058699	0.37597	1.043	2.8297	0.44625 ];

particlesB6 = [ -0.095913	-0.1817	-0.10521	0.43851	0.20661	-1	0.58426	0.61941	1.1791	2.9029	0.41225;
    -0.058671	-0.35286	0.49664	0.12088	0.31666	-1	-0.13826	-0.43334	1.1747	2.9165	0.41737;
    -0.10565	0.34853	-0.54094	0.25548	0.12673	0.52359	-1	1	1.1838	3.0477	0.41451 ];

particlesB8 = [ -0.55738	0.6219	-0.31968	0.045746	-0.070674	0.66418	-0.35341	-0.45893	1.2848	2.984	0.37425;
    -0.60114	0.4061	0.38183	0.044342	-0.005801	0.69928	-0.77527	-0.78553	1.2649	2.9687	0.38765;
    -0.16974	0.12319	0.57286	-0.47868	-0.37655	0.34343	-0.22135	0.059018	1.2496	2.6008	0.40282 ];

particlesB10 = [ -0.044198	-0.015433	-0.14805	0.32411	-0.088617	0.011836	1	-0.7328	1.2807	2.8609	0.36592;
    0.53415	-0.46122	-0.95307	0.76981	-0.22034	0.93779	-0.60791	-0.32023	1.3263	2.8439	0.35271;
    -0.07758	0.032777	0.39127	-0.082905	-0.32245	0.61233	-0.10411	-0.32031	1.2825	2.6859	0.37292 ];






JBestValuesB2 = [1.0837 , 1.084 1.0834 ];
JBestValuesB4 = [ 1.486,1.488, 1.4893 ];
JBestValuesB6 = [ 1.5914 1.5921 1.5983 ];
JBestValuesB8 = [ 1.6591 1.6526 1.6525 ];
JBestValuesB10 = [ 1.6466 1.6791 1.6554 ];




particles = particlesB2;
JBestValues = JBestValuesB2;

partSize = size(particles);
numParticles = partSize(1);

% [ xi0, xi1, xi2, xi3, v0, v1, v2, v3, v4, t1, e, t2 ]  

for numPart = 1:numParticles
    
    deltaT1 = particles(numPart, 9);
    t1 = 0:.01:deltaT1;
    deltaT2 = particles(numPart, 11);
    tminust2 = 0:.01:deltaT2;
    
    xi0 = particles(numPart, 1);
    xi1 = particles(numPart, 2);
    xi2 = particles(numPart, 3);
    xi3 = particles(numPart, 4);
    
    v0 = particles(numPart, 5);
    v1 = particles(numPart, 6);
    v2 = particles(numPart, 7);
    v3 = particles(numPart, 8);
    
    delta1 = xi0+xi1.*t1+xi2.*t1.^2+xi3.*t1.^3;
    delta2 = v0+v1.*tminust2+v2.*tminust2.^2+v3.*tminust2.^3;
    
    subplot(2,1,1);
    plot(t1, rad2deg(delta1), '-o');
    title("Thrust arc 1 angle B=10");
    ylabel("Angle (radians)");
    xlabel("Time (TU)");
        legend(['JBest=' num2str(JBestValues(1))],['JBest=' num2str(JBestValues(2))], ['JBest=' num2str(JBestValues(3))]);


    hold on;
    
    subplot(2,1,2);
    plot(tminust2, rad2deg(delta2), '-o');
    title("Thrust arc 2 angle B=4");
    ylabel("Angle (radians)");
    xlabel("Time (TU)");
        legend(['JBest=' num2str(JBestValues(1))],['JBest=' num2str(JBestValues(2))], ['JBest=' num2str(JBestValues(3))]);


    hold on;
    
    
    
    
end



