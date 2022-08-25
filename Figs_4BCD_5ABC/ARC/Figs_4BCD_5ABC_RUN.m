%% Code for Figs 4BCD,5ABC: the optimal number of sentinels to include in the sample and the corresponding reductions in EDP
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 17th February 2022

%% -----------------------------------------------------------------------------------------------
% For a given number of sentinels added to the population, this code iterates over
% specified values of the sample size and sample interval and in each instance applies a
% Bayesian optimisation algorithm to compute the optimal number of sentinels to include in
% the sample and the corresponding reduction in EDP compared to the baseline level.

function Figs_4BCD_5ABC_RUN(ID)

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
sampleSizeVec = 25:25:200;
% Vector of sample intervals to consider
sampleIntervalVec = 30:30:150;

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
savePath = './Figs_3BCD_4ABC_results/';
% Define random number generator
rng('shuffle');

%% ------------------------------------------------------------------------
% END USER INPUT

%% ------------------------------------------------------------------------
% PRELIMINARY CALCULATIONS

% Total population size (P)
P = Pc+Ps;

% Extract sample size and sample interval from lists
[p,q] = meshgrid(sampleSizeVec,sampleIntervalVec);
pairs = [p(:) q(:)];
pairs_choose = pairs(ID,:);
sampleSize = pairs_choose(1);
sampleInterval = pairs_choose(2);

% Compute maximum allowable number of sentinels in sample
maxNoSentinels = min(Ps,sampleSize);

%% ------------------------------------------------------------------------
% COMPUTE BASELINE EDP
simData = runSpread_1(Pc,D0,U0,betaC,epsilonC,gammaC,numSims,tFinal,progress);
[~, EDP, ~, ~, ~] = runSampling_1(simData,numRuns,Pc,sampleSize,sampleInterval,tFinal,progress);
baselineEDP = EDP;

%% ------------------------------------------------------------------------
% PERFORM BAYESIAN OPTIMISATION

if maxNoSentinels == 0 % then no sentinels allowed, so optimum is 0
    
    results = [sampleSize, sampleInterval, 0, 0, baselineEDP];
    
else % Run spread simulations in the presence of sentinels

    simData = runSpread_2(P,Ps,D0,U0,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,numSims,tFinal,progress);
    % Define the number of sentinels in the sample as an optimisable variable
    numSentinelsSampled = optimizableVariable('S',[0,maxNoSentinels],'Type','integer');

    % Specify objective function input
    fun = @(z)objFunc(simData,numRuns,P,Ps,sampleSize-z.S,z.S,sampleInterval,tFinal,progress,baselineEDP);
    % Perform Bayesian optimisation
    optResults = bayesopt(fun,numSentinelsSampled,'IsObjectiveDeterministic',false,'MaxObjectiveEvaluations',numIterations,'AcquisitionFunctionName','expected-improvement-plus','PlotFcn',[]);

    % Extract optimal parameters
    optimalParams = optResults.XAtMinEstimatedObjective;
    % Extract optimal number of sentinels to include in sample
    optimalSent = optimalParams.S;

    % Extract value of objective function at optimum
    valueAtOpt = optResults.MinEstimatedObjective;
    % Compute resultant EDP at optimum
    resultantEDP = 0.01*valueAtOpt*baselineEDP+baselineEDP;

    % Store results
    results = [sampleSize, sampleInterval, optimalSent, valueAtOpt, resultantEDP];
    
end

%% ------------------------------------------------------------------------
% WRITE RESULTS TO .TXT FILES

% Make new directory to store results
mkdir(savePath)

T = table(Pc,Ps,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,D0,U0,numSims,numRuns,tFinal,numIterations);
if ID==1
    writetable(T,[savePath 'params.txt'],'Delimiter','tab');
end

filename = sprintf('results_%d.txt',ID);
writematrix(results,[savePath filename],'Delimiter','tab');

end
