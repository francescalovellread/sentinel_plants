%% Code for Figs 4BCD,5ABC: the optimal number of sentinels to include in the sample and the corresponding reductions in EDP
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 4th August 2022

%% -----------------------------------------------------------------------------------------------
% For a given number of sentinels added to the population, this code iterates over
% specified values of the sample size and sample interval and in each instance applies a
% Bayesian optimisation algorithm to compute the optimal number of sentinels to include in
% the sample and the corresponding reduction in EDP compared to the baseline level.
clear; close all;

%% ------------------------------------------------------------------------
% DEFINE MODEL PARAMETERS

% Crop population size (PC)
Pc = 1000; 
% Sentinel population size (PS)
Ps = 0; 
% Transmission coefficient for 'Detectable' crops (betaC) 
r=0.05; betaC = r/Pc;
% Transmission coefficient for 'Detectable' sentinels (betaS)
betaS = betaC;
% Transmission scaling factor for 'Undetectable' crops (epsilonC)
epsilonC = 0.015; 
% Transmission scaling factor for 'Undetectable' sentinels (epsilonS)
epsilonS = 0.1;
% Duration of crop 'Undetectable' period (gammaC)
gammaC = 452;
% Duration of sentinel 'Undetectable' period (gammaS)
gammaS = 49;
% Initial numbers of 'Detectable' and 'Undetectable' plants
D0 = 0; U0 = 1; 

% Vector of sample sizes to consider
sampleSizeVec = 20:20:100;
% Vector of sample intervals to consider
sampleIntervalVec = 30:30:120;

% Number of spread simulations to perform
numSims = 100; 
% Number of sampling simulations to perform
numRuns = 100;
% Maximum run time for spread and sampling simulations
tFinal = 5000; 
% Specify whether to display progress messages
progress = "no";
% Number of iterations for Bayesian optimisation algorithm
numIterations = 30;

%% ------------------------------------------------------------------------
% END USER INPUT

%% ------------------------------------------------------------------------
% PRELIMINARY CALCULATIONS

% Total population size (P)
P = Pc+Ps;
% Extract lengths of sampleSizeVec and sampleIntervalVec
nSize = length(sampleSizeVec); 
nInterval = length(sampleIntervalVec);

%% ------------------------------------------------------------------------
% ITERATE OVER ALL SAMPLE SIZES AND SAMPLE INTERVALS AND PERFORM BAYESIAN OPTIMISATION

% Create empty matrices to store generated values
optimalSentStore = zeros(nInterval,nSize); % for the optimum number of sentinels in sample
valueAtOptStore = zeros(nInterval,nSize); % for the corresponding percentage change in EDP at the optimum
resultantEDPStore = zeros(nInterval,nSize); % for the EDP at the optimum

for j = 1:nSize % Iterate over sample sizes
    sampleSize = sampleSizeVec(j);
    
    for k = 1:nInterval % Iterate over sample intervals
        sampleInterval = sampleIntervalVec(k);

        % Compute baseline EDP
        simData = runSpread_1(Pc,D0,U0,betaC,epsilonC,gammaC,numSims,tFinal,progress);
        [sampleData, EDP, EDT, EDR, perc95] = runSampling_1(simData,numRuns,Pc,sampleSize,sampleInterval,tFinal,progress);
        baselineEDP = EDP;

        % Compute maximum allowable number of sentinels in sample
        maxNoSentinels = min(Ps,sampleSize);

        if maxNoSentinels == 0 % then no sentinels allowed, so optimum is 0
            optimalSentStore(k,j) = 0;
            valueAtOptStore(k,j) = 0;
            resultantEDPStore(k,j) = baselineEDP;
            
        else
            % Run spread simulations in the presence of sentinels
            simData = runSpread_2(P,Ps,D0,U0,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,numSims,tFinal,progress);
            % Define the number of sentinels in the sample as an optimisable variable
            numSentinelsSampled = optimizableVariable('S',[0,maxNoSentinels],'Type','integer');
            % Specify objective function input (function defined at end of script)
            fun = @(z)objFunc(simData,numRuns,P,Ps,sampleSize-z.S,z.S,sampleInterval,tFinal,progress,baselineEDP);
            % Perform Bayesian optimisation
            optResults = bayesopt(fun,numSentinelsSampled,'IsObjectiveDeterministic',false,'MaxObjectiveEvaluations',numIterations,'AcquisitionFunctionName','expected-improvement-plus','PlotFcn',[]);

            % Extract optimal number of sentinels to include in sample
            optimalParams = optResults.XAtMinEstimatedObjective;
            optimalSent = optimalParams.S;
            % Extract value of objective function at optimum
            valueAtOpt = optResults.MinEstimatedObjective;
            % Compute resultant EDP at optimum
            resultantEDP = 0.01*valueAtOpt*baselineEDP+baselineEDP;
            % Store results
            optimalSentStore(k,j) = optimalSent;
            valueAtOptStore(k,j) = valueAtOpt;
            resultantEDPStore(k,j) = resultantEDP;
            
        end
    end
end

%% ------------------------------------------------------------------------
% PLOT

% Define x plot values
xplot = sampleSizeVec;

% Define colours for IBM colourblind safe palette
IBMblue = [100,143,255]/256;
IBMpurple = [120,93,241]/256;
IBMpink = [221,37,128]/256;
IBMorange = [254,97,0]/256;
IBMyellow = [255,176,0]/256;

%% ------------------------------------------------------------------------
% PLOT OPTIMAL NUMBER OF SENTINELS TO INCLUDE IN SAMPLE

figure(); grid on; hold on;

% Plot line marking upper bound on the number of sentinels sampled
boundline1 = line([sampleSizeVec(1),Ps],[sampleSizeVec(1),Ps]);
boundline2 = line([Ps,sampleSizeVec(end)],[Ps,Ps]);
boundline1.LineWidth = 2; boundline2.LineWidth = 2;
boundline1.Color = [0 0 0]; boundline2.Color = [0 0 0];

% Plot optimal nuber of sentinels to include over range of sample sizes for each sample
% interval considered
myplot(1) = plot(xplot,optimalSentStore(1,:));
myplot(2) = plot(xplot,optimalSentStore(2,:));
myplot(3) = plot(xplot,optimalSentStore(3,:));
myplot(4) = plot(xplot,optimalSentStore(4,:));

% Title
title(['P_S =' 32 num2str(Ps) 32 'sentinels']);

% Format x axis
xlabel('Sample size (N)');
xlim([sampleSizeVec(1) sampleSizeVec(end)]);

% Format y axis
ylabel('Optimal number of sentinels in sample (N_S^*)');
ylim([0 100]);
yticks = 0:10:100;
set(gca, 'YTick', yticks, 'Fontsize', 16);

% Legend
legEntries = {'min(P_S,N)','\Delta = 30 days','\Delta = 60 days','\Delta = 90 days','\Delta = 120 days'};
leg = legend([boundline1 myplot(1) myplot(2) myplot(3) myplot(4)],legEntries);
leg.Location = 'northwest';

% Define plot markers, colours and styles
myplot(1).Marker = '^';
myplot(2).Marker = '*'; 
myplot(3).Marker = 'o'; 
myplot(4).Marker = 'x'; 

myplot(1).Color = IBMblue;
myplot(1).MarkerFaceColor = IBMblue;
myplot(1).MarkerEdgeColor = IBMblue;

myplot(2).Color = IBMpink;
myplot(2).MarkerFaceColor = IBMpink;
myplot(2).MarkerEdgeColor = IBMpink;

myplot(3).Color = IBMorange;
myplot(3).MarkerFaceColor = IBMorange;
myplot(3).MarkerEdgeColor = IBMorange;

myplot(4).Color = IBMpurple;
myplot(4).MarkerFaceColor = IBMpurple;
myplot(4).MarkerEdgeColor = IBMpurple;

for l=1:4;  myplot(l).LineWidth = 2; myplot(l).LineStyle = ':'; myplot(l).MarkerSize = 10; end

box off; set(gca,'Fontsize',16,'Linewidth',2);

%% ------------------------------------------------------------------------
% PLOT CORRESPONDING REDUCTIONS IN EDP

figure(); grid on; hold on;

xplot = sampleSizeVec;

% Plot EDP reductions at optimum
myplot(1) = plot(xplot,valueAtOptStore(1,:));
myplot(2) = plot(xplot,valueAtOptStore(2,:));
myplot(3) = plot(xplot,valueAtOptStore(3,:));
myplot(4) = plot(xplot,valueAtOptStore(4,:));

% Title
title('P_S = 50 sentinels');

% Format x axis
xlabel('Sample size (N)');
xlim([sampleSizeVec(1) sampleSizeVec(end)]); 

% Format y axis
ylabel('Change in EDP from baseline at optimum (%)');

% Legend
legEntries = {'\Delta = 30 days','\Delta = 60 days','\Delta = 90 days','\Delta = 120 days'};
legend([myplot(1) myplot(2) myplot(3) myplot(4)],legEntries);

% Define plot markers, colours and styles
myplot(1).Marker = '^';
myplot(2).Marker = '*'; 
myplot(3).Marker = 'o'; 
myplot(4).Marker = 'x'; 

myplot(1).Color = IBMblue;
myplot(1).MarkerFaceColor = IBMblue;
myplot(1).MarkerEdgeColor = IBMblue;

myplot(2).Color = IBMpink;
myplot(2).MarkerFaceColor = IBMpink;
myplot(2).MarkerEdgeColor = IBMpink;

myplot(3).Color = IBMorange;
myplot(3).MarkerFaceColor = IBMorange;
myplot(3).MarkerEdgeColor = IBMorange;

myplot(4).Color = IBMpurple;
myplot(4).MarkerFaceColor = IBMpurple;
myplot(4).MarkerEdgeColor = IBMpurple;

for l=1:4;  myplot(l).LineWidth = 2; myplot(l).LineStyle = ':'; myplot(l).MarkerSize = 10; end

box off; set(gca,'Fontsize',16,'Linewidth',2);

%% ------------------------------------------------------------------------
% WRITE RESULTS TO TXT FILES
% dlmwrite('sampleSizeVec.txt',sampleSizeVec)
% dlmwrite('sampleIntervalVec.txt',sampleIntervalVec)
% dlmwrite(strcat('valueAtOptStore_',num2str(numSentinels),'.txt'),valueAtOptMat)

%% ------------------------------------------------------------------------
% DEFINE OBJECTIVE FUNCTION TO BE MINIMISED BY BAYESOPT

function OUTPUT = objFunc(simData,numRuns,popSize,numSentinels,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,~,baselineEDP)
    [~, ~, ~, ~, ECDP, ~] = runSampling_2(simData,numRuns,popSize,numSentinels,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,"no");
    OUTPUT = 100*(ECDP-baselineEDP)/baselineEDP;
end