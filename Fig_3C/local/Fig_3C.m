%% Code for Fig 3C: computing the baseline EDP for a range of sample sizes and intervals
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 4th August 2022

%% -----------------------------------------------------------------------------------------------
% This code uses the functions 'runSpread_1' and 'runSampling_1' to compute the
% baseline EDP across a specified range of sample sizes and intervals. 
clear; close all;

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
sampleSizeVec = 25:25:200;
% Vector of sample intervals to consider
sampleIntervalVec = 30:30:150;

% Number of spread simulations to perform
numSims = 20000; 
% Number of sampling simulations to perform
numRuns = 20000;
% Maximum run time for spread and sampling simulations
tFinal = 5000; 
% Specify whether to display progress messages
progress = "no";

%% ------------------------------------------------------------------------
% END USER INPUT

%% ------------------------------------------------------------------------
% GENERATE EPIDEMIC CURVES

simData = runSpread_1(P,D0,U0,beta,epsilon,gamma,numSims,tFinal,progress);

%% ------------------------------------------------------------------------
% PERFORM SAMPLING SCHEME FOR RANGE OF SAMPLE SIZES AND INTERVALS

% Create empty matrix to store generated values
EDPstore = zeros(length(sampleSizeVec),length(sampleIntervalVec));

for i=1:length(sampleSizeVec) % Iterate over sample sizes
    sampleSize = sampleSizeVec(i);
    
    for j=1:length(sampleIntervalVec) % Iterate over sample intervals
        sampleInterval = sampleIntervalVec(j);
        
        % Print progress message if specified
        if progress == "yes"; fprintf(strcat('Sample size =',32,num2str(sampleSize),',',32,'Sample interval =',32,num2str(sampleInterval),'\n')); end

        % Perform sampling simulations and store results
        [sampleData, EDP, EDT, EDR, perc95] = runSampling_1(simData,numRuns,P,sampleSize,sampleInterval,tFinal,progress);
        EDPstore(i,j) = EDP;     
    end
end

%% ------------------------------------------------------------------------
% PLOT

figure();

% Create contour plot of EDPs expressed as percentage of total population
[C,h] = contourf(100*EDPstore'/P);

% Set contour and label styles and specify manual labelling
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

% Define colourbar
colbar = colorbar;
ylabel(colbar, 'Baseline EDP (% of population)', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

% Set x axis label and ticks
xlabel('Sample size (N=N_C)')
xticks = 1:1:length(sampleSizeVec);
xticklabels = strsplit(num2str(sampleSizeVec));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);

% Set y axis label and ticks
ylabel('Sample interval (\Delta days)')
yticks = 1:1:length(sampleIntervalVec);
yticklabels = strsplit(num2str(sampleIntervalVec));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);

box off; set(gca,'Fontsize',16,'Linewidth',2);
