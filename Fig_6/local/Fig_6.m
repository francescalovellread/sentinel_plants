%% Code for Fig 6: Optimising the total number of sentinels in the population and the sample
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 4th August 2022

%% -----------------------------------------------------------------------------------------------
% This code iterates over specified values of the sample size and sample interval,
% and in each instance applies a Bayesian optimisation algorithm to compute the
% optimal number of sentinels to include in the population and in the sample. It
% also computes the corresponding reduction in EDP compared to the baseline level.
clear; close all;

%% ------------------------------------------------------------------------
% DEFINE MODEL PARAMETERS

% Crop population size (PC)
Pc = 1000; 
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
sampleSizeVec = 100:50:200;
% Vector of sample intervals to consider
sampleIntervalVec = 30:30:150;

% Number of spread simulations to perform
numSims = 1000; 
% Number of sampling simulations to perform
numRuns = 1000;
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

nSize = length(sampleSizeVec); 
nInterval = length(sampleIntervalVec);

%% ------------------------------------------------------------------------
% ITERATE OVER ALL SAMPLE SIZES AND SAMPLE INTERVALS AND PERFORM BAYESIAN OPTIMISATION

% Create empty matrix to store results
resultsStoreCell = cell(nSize,nInterval);

for i=1:nSize % Iterate over sample sizes
    sampleSize = sampleSizeVec(i);
    
    for j=1:nInterval % Iterate over sample intervals
        sampleInterval = sampleIntervalVec(j);
        
        % Compute baseline EDP
        simData = runSpread_1(Pc,D0,U0,betaC,epsilonC,gammaC,numSims,tFinal,progress);
        [sampleData, EDP, EDT, EDR, perc95] = runSampling_1(simData,numRuns,Pc,sampleSize,sampleInterval,tFinal,progress);
        baselineEDP = EDP;

        % Define the number of sentinels in the population as an optimisable variable
        totalNumSentinels = optimizableVariable('Stot',[0,350],'Type','integer');
        % Define the number of sentinels in the sample as an optimisable variable        
        numSentinelsSampled = optimizableVariable('Ssamp',[0,sampleSize],'Type','integer');

        % Specify objective function input (function defined at end of script)
        fun = @(z)objFunc(Pc+z.Stot,z.Stot,D0,U0,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,numSims,tFinal,progress,numRuns,sampleSize-z.Ssamp,z.Ssamp,sampleInterval,baselineEDP);
        % Perform Bayesian optimisation
        optResults = bayesopt(fun,[totalNumSentinels,numSentinelsSampled],'IsObjectiveDeterministic',false,'MaxObjectiveEvaluations',numIterations,'AcquisitionFunctionName','expected-improvement-plus','XConstraintFcn',@xconstraint,'UseParallel',false,'PlotFcn',[]);
       
        % Extract optimal parameters
        optimalParams = optResults.XAtMinEstimatedObjective;
        % Extract optimal number of sentinels to include in population
        optSentinelsAdded = optimalParams.Stot;
        % Extract optimal number of sentinels to include in sample
        optSentinelsSampled = optimalParams.Ssamp;
        
        % Extract value of objective function at optimum
        valueAtOpt = optResults.MinEstimatedObjective;
        
        % Compute resultant EDP at optimum
        resultantEDP = 0.01*valueAtOpt*baselineEDP+baselineEDP;

        % Store results
        resultsVec = [sampleSize, sampleInterval, optSentinelsAdded, optSentinelsSampled, valueAtOpt, baselineEDP, resultantEDP];
        resultsStoreCell{i,j} = resultsVec;
    end
end

% Reformat stored results into matrix
resultsStore = [];
for i=1:nSize
    for j=1:nInterval
        resultsStore = [resultsStore; resultsStoreCell{i,j}];
    end
end

% Extract matrix of the optimal number of sentinels added to the population
optTotSentinels = reshape(resultsStore(:,3),[nInterval nSize]);

% Extract matrix of the optimal proportion of sentinels in the sample
optSentinelsSampled = reshape(resultsStore(:,4),[nInterval nSize]);
optSentProp = optSentinelsSampled./min(optTotSentinels,sampleSizeVec);

% Extract matrix of the reduction in EDP at the optimum
EDPreduction = reshape(resultsStore(:,5),[nInterval nSize]);

% Extract matrix of the resultant EDP at the optimum
EDPresultant = reshape(resultsStore(:,7),[nInterval nSize]);
EDPresultant = 100*EDPresultant/Pc;

% Compute matrix of the relative sentinel utility
utility = -EDPreduction./optTotSentinels;

%% ------------------------------------------------------------------------
% PLOT OPTIMAL NUMBER OF SENTINELS

figure();

% Create contour plot of optimal number of sentinels
[C,h] = contourf(optTotSentinels);

% Set contour and label styles and specify manual labelling
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

% Define colourbar
colbar = colorbar;
ylabel(colbar, 'Optimal number of sentinels (P_S^*)', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

% Format x axis
xlabel('Sample size (N)');
xticks = 1:1:5;
xticklabels = strsplit(num2str(sampleSizeVec));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

% Format y axis
ylabel('Sample interval (\Delta days)');
yticks = 1:1:4;
yticklabels = strsplit(num2str(sampleIntervalVec));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

set(gca,'fontsize',16,'Linewidth',2); box off;

%% ------------------------------------------------------------------------
% PLOT OPTIMAL SENTINEL PROPORTION

figure(); 

% Create contour plot of optimal sentinel proportion
[C,h] = contourf(optSentProp);

% Set contour and label styles and specify manual labelling
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end
caxis([0 1]);

% Define colourbar
colbar = colorbar;
ylabel(colbar, 'Optimal proportion of sentinels in sample', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

% Format x axis
xlabel('Sample size (N)');
xticks = 1:1:5;
xticklabels = strsplit(num2str(sampleSizeVec));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

% Format y axis
ylabel('Sample interval (\Delta days)');
yticks = 1:1:4;
yticklabels = strsplit(num2str(sampleIntervalVec));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

set(gca,'fontsize',16,'Linewidth',2); box off;

%% ------------------------------------------------------------------------
% PLOT REDUCTION IN EDP

figure(); 

% Create contour plot of EDP reduction
[C,h] = contourf(EDPreduction);

% Set contour and label styles and specify manual labelling
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

% Define colourbar
colbar = colorbar;
ylabel(colbar, 'Change in EDP from baseline at optimum (%)', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

% Format x axis
xlabel('Sample size (N)');
xticks = 1:1:5;
xticklabels = strsplit(num2str(sampleSizeVec));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

% Format y axis
ylabel('Sample interval (\Delta days)');
yticks = 1:1:4;
yticklabels = strsplit(num2str(sampleIntervalVec));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

set(gca,'fontsize',16,'Linewidth',2); box off;

%% ------------------------------------------------------------------------
% PLOT RESULTANT EDP

figure();

% Create contour plot of resultant EDP
[C,h] = contourf(EDPresultant);

% Set contour and label styles and specify manual labelling
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

% Define colourbar
colbar = colorbar;
ylabel(colbar, 'EDP at optimum (% of crop population)', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

% Format x axis
xlabel('Sample size (N)');
xticks = 1:1:5;
xticklabels = strsplit(num2str(sampleSizeVec));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

% Format y axis
ylabel('Sample interval (\Delta days)');
yticks = 1:1:4;
yticklabels = strsplit(num2str(sampleIntervalVec));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

set(gca,'fontsize',16,'Linewidth',2); box off;

%% ------------------------------------------------------------------------
% PLOT SENTINEL UTILITY

figure(); 

% Create contour plot of sentinel utility
[C,h] = contourf(utility);

% Set contour and label styles and specify manual labelling
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

% Define colourbar
colbar = colorbar;
ylabel(colbar, 'Relative utility of sentinel', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

% Format x axis
xlabel('Sample size (N)');
xticks = 1:1:nSize;
xticklabels = strsplit(num2str(sampleSizeVec));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

% Format y axis
ylabel('Sample interval (\Delta days)');
yticks = 1:1:nInterval;
yticklabels = strsplit(num2str(sampleIntervalVec));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

set(gca,'fontsize',16,'Linewidth',2); box off;

%% ------------------------------------------------------------------------
% DEFINE OBJECTIVE FUNCTION TO BE MINIMISED BY BAYESOPT

function OUTPUT = objFunc(popSize,numSentinels,initialI,initialC,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,numSims,tFinal,progress,numRuns,cropSampleSize,sentinelSampleSize,sampleInterval,baselineEDP)
    if numSentinels==0
        OUTPUT = 0;
    else
        simData = runSpread_2(popSize,numSentinels,initialI,initialC,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,numSims,tFinal,progress);
        [~, ~, ~, ~, ECDP, ~] = runSampling_2(simData,numRuns,popSize,numSentinels,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,progress); 
        OUTPUT = 100*(ECDP-baselineEDP)/baselineEDP;
    end
end

%% ------------------------------------------------------------------------
% DEFINE DETERMINISTIC CONSTRAINT FUNCTION FOR BAYESOPT

function tf = xconstraint(X)
    tf = X.Ssamp<=X.Stot;
end
