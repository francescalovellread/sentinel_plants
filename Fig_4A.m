%% Code for Fig 4A: Example of computing the optimal value of N_S
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 4th August 2022

%% -----------------------------------------------------------------------------------------------
% This code demonstrates finding the optimal number of sentinels to include in the sample
% for given values of the number of sentinels added (P_S), the sample size (N) and the
% sample interval (Delta), using a Bayesian optimisation algorithm.
clear; close all;

%% ------------------------------------------------------------------------
% DEFINE MODEL PARAMETERS

% Crop population size (PC)
Pc = 1000; 
% Sentinel population size (PS)
Ps = 50; 
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
% Sample size
sampleSize = 50;
% Sample interval
sampleInterval = 30;
% Number of spread simulations to perform
numSims = 10000; 
% Number of sampling simulations to perform
numRuns = 10000;
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
% Maximum possible number of sentinels in sample
maxNoSentinels = min(Ps,sampleSize);

%% ------------------------------------------------------------------------
% COMPUTE THE BASELINE EDP IN THE ABSENCE OF SENTINELS

simData = runSpread_1(Pc,D0,U0,betaC,epsilonC,gammaC,numSims,tFinal,progress);
[sampleData, EDP, EDT, EDR, perc95] = runSampling_1(simData,numRuns,Pc,sampleSize,sampleInterval,tFinal,progress);
baselineEDP = EDP;

%% ------------------------------------------------------------------------
% RUN SPREAD SIMULATIONS IN THE PRESENCE OF SENTINELS

simData = runSpread_2(P,Ps,D0,U0,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,numSims,tFinal,progress);

%% ------------------------------------------------------------------------
% PERFORM BAYESIAN OPTMISATION ON THE SENTINEL SAMPLE SIZE Ns

% Define the number of sentinels in the sample as an optimisable variable
numSentinelsSampled = optimizableVariable('Ns',[0,maxNoSentinels],'Type','integer');
% Specify objective function input (function defined at end of script)
fun = @(z)objFunc(simData,numRuns,P,Ps,sampleSize-z.Ns,z.Ns,sampleInterval,tFinal,progress,baselineEDP);
% Perform Bayesian optimisation
optResults = bayesopt(fun,numSentinelsSampled,'IsObjectiveDeterministic',false,'MaxObjectiveEvaluations',numIterations,'AcquisitionFunctionName','expected-improvement-plus');

%% ------------------------------------------------------------------------
% EXTRACT RESULTS FROM BAYESIAN OPTIMISATION

% Extract optimal number of sentinels to include in sample
optimalParams = optResults.XAtMinEstimatedObjective;
optimalSent = optimalParams.Ns;
% Extract value of objective function at optimum
valueAtOpt = optResults.MinEstimatedObjective;
% Extract model mean
figure(1); 
kids = get(gca,'Children');
modelmean = kids(7);
yvals = unique(modelmean.YData,'stable');
xvals = 0:maxNoSentinels;

%% ------------------------------------------------------------------------
% PLOT

figure(); hold on; grid on;

% Plot model mean
meanplot = plot(xvals,yvals);
plotcol = [108 197 20]/256;
meanplot.Color = plotcol; meanplot.LineWidth = 2;
meanplot.LineStyle = ':'; meanplot.Marker = 'x';

% Plot minimum point
minPoint = plot(optimalSent,valueAtOpt); 
minPoint.Marker = 'o'; minPoint.MarkerSize = 12;
minPoint.MarkerFaceColor = plotcol; minPoint.MarkerEdgeColor = plotcol;
box off; set(gca,'Fontsize',16,'Linewidth',2);

% Plot baseline
bl = yline(0);
bl.LineStyle = '--'; bl.LineWidth = 2; bl.Color = [0 0 0];

title(strcat('P_S =',32,num2str(Ps),',',32,'N =',32,num2str(sampleSize),',',32,'\Delta =',32,num2str(sampleInterval),32,'days'));
xlabel('Number of sentinels in sample (N_S)')
ylabel('Change in EDP from baseline (%)')

%% ------------------------------------------------------------------------
% DEFINE OBJECTIVE FUNCTION TO BE MINIMISED BY BAYESOPT

function OUTPUT = objFunc(simData,numRuns,popSize,numSentinels,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,~,baselineEDP)
    [~, ~, ~, ~, ECDP, ~] = runSampling_2(simData,numRuns,popSize,numSentinels,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,"no");
    OUTPUT = 100*(ECDP-baselineEDP)/baselineEDP;
end