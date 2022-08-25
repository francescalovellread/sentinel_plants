%% Code for Fig 3C: computing the baseline EDP for a range of sample sizes and intervals
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 4th August 2022

%% -----------------------------------------------------------------------------------------------
% This code uses the functions 'runSpread_1' and 'runSampling_1' to compute the
% baseline EDP across a specified range of sample sizes and intervals. 

function Fig_3C_RUN(ID)

%% ------------------------------------------------------------------------
% DEFINE MODEL PARAMETERS

% Population size (P=P_C)
P = 1000; 
% Transmission coefficient for 'Detectable' crops (beta_C)
r=0.05; beta = r/P;
% Transmission scaling factor for 'Undetectable' crops (epsilon_C)
epsilon = 0.015; 
% Duration of crop 'Undetectable' period (gamma_C)
gamma = 452;
% Initial numbers of 'Detectable' and 'Undetectable' crops
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
% Define file path for save location
savePath = './Fig_3C_results/';
% Define random number generator
rng('shuffle');

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
% GENERATE EPIDEMIC CURVES

simData = runSpread_1(P,D0,U0,beta,epsilon,gamma,numSims,tFinal,progress);

%% ------------------------------------------------------------------------
% PERFORM SAMPLING SCHEME FOR RANGE OF SAMPLE SIZES AND INTERVALS

[~, EDP, ~, ~, ~] = runSampling_1(simData,numRuns,P,sampleSize,sampleInterval,tFinal,progress);    
results = [sampleSize, sampleInterval, EDP];
    
%% ------------------------------------------------------------------------
% WRITE RESULTS TO .TXT FILES

% Make new directory to store results
mkdir(savePath)

T = table(P,r,beta,epsilon,gamma,D0,U0,numSims,numRuns,tFinal);
if ID==1
    writetable(T,[savePath 'params.txt'],'Delimiter','tab');
end

filename = sprintf('results_%d.txt',ID);
writematrix(results,[savePath filename],'Delimiter','tab');

end