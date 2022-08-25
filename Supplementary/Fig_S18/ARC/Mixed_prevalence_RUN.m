%% Code for Supplementary Fig S18: Considering a mixture of crop and sentinel detection prevalences
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 4th August 2022

%% -----------------------------------------------------------------------------------------------
% This code is equivalent to the code for Fig 6 of the main text, but considers the case
% in which the detection prevalence amongst sentinel plants is also accounted for (i.e.
% the objective function is EDP + (omega_S*EDP_sent).

function Mixed_prevalence_RUN(ID)

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
sampleSizeVec = 25:5:200;
% Vector of sample intervals to consider
sampleIntervalVec = 30:5:150;

% Number of spread simulations to perform
numSims = 25000; 
% Number of sampling simulations to perform
numRuns = 25000;
% Maximum run time for spread and sampling simulations
tFinal = 5000; 
% Specify whether to display progress messages
progress = "no";
% Number of iterations for Bayesian optimisation algorithm
numIterations = 30;
% Specify number of parallel workers
nWorkers = 1;
% Define file path for save location
savePath = './Fig_5_results/';
% Define random number generator
rng('shuffle');

% Specify relative weight of sentinel detection prevalence
omegaS = 0.5;

%% ------------------------------------------------------------------------
% END USER INPUT

%% ------------------------------------------------------------------------
% EXTRACT SAMPLE SIZE AND SAMPLE INTERVAL FROM LISTS

[p,q] = meshgrid(sampleSizeVec,sampleIntervalVec);
pairs = [p(:) q(:)];
pairs_choose = pairs(ID,:);
sampleSize = pairs_choose(1);
sampleInterval = pairs_choose(2);

%% ------------------------------------------------------------------------
% COMPUTE BASELINE EDP
simData = runSpread_1(Pc,D0,U0,betaC,epsilonC,gammaC,numSims,tFinal,progress);
[~, EDP, ~, ~, ~] = runSampling_1(simData,numRuns,Pc,sampleSize,sampleInterval,tFinal,progress);
baselineEDP = EDP;

%% ------------------------------------------------------------------------
% PERFORM BAYESIAN OPTIMISATION

% Define the number of sentinels in the population as an optimisable variable
totalNumSentinels = optimizableVariable('Stot',[0,350],'Type','integer');
% Define the number of sentinels in the sample as an optimisable variable        
numSentinelsSampled = optimizableVariable('Ssamp',[0,sampleSize],'Type','integer');

% Specify objective function input 
fun = @(z)objFunc(Pc+z.Stot,z.Stot,D0,U0,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,numSims,tFinal,progress,numRuns,sampleSize-z.Ssamp,z.Ssamp,sampleInterval,baselineEDP,omegaS);
% Perform Bayesian optimisation
optResults = bayesopt(fun,[totalNumSentinels,numSentinelsSampled],'IsObjectiveDeterministic',false,'MaxObjectiveEvaluations',numIterations,'AcquisitionFunctionName','expected-improvement-plus','XConstraintFcn',@objConstraint,'UseParallel',false,'PlotFcn',[]);

%% ------------------------------------------------------------------------
% EXTRACT AND STORE RESULTS

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
results = [sampleSize, sampleInterval, optSentinelsAdded, optSentinelsSampled, valueAtOpt, baselineEDP, resultantEDP];

%% ------------------------------------------------------------------------
% WRITE RESULTS TO .TXT FILES

% Make new directory to store results
mkdir(savePath)

T = table(Pc,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,D0,U0,numSims,numRuns,tFinal,numIterations);
if ID==1
    writetable(T,[savePath 'params.txt'],'Delimiter','tab');
end

filename = sprintf('results_%d.txt',ID);
writematrix(results,[savePath filename],'Delimiter','tab');

end