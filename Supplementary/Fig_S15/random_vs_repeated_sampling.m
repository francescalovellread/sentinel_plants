%% Code for Supplementary Fig S15: Demonstrating the difference between random and repeated sampling
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 4th August 2022

%% -----------------------------------------------------------------------------------------------
% This code performs spread and sampling simulations on a simple SI model for a range of
% sample sizes and sample intervals. It performs both random sampling and repeated sampling
% to demonstrate how repeated sampling can lead to a worse outcome than random sampling (generates
% Supplementary Fig 15).
clear; close all;

%% ------------------------------------------------------------------------
% DEFINE MODEL PARAMETERS

% Population size
popSize = 1000; 
% Transmission coefficient
beta=5e-6; 
% Initial number of infected individuals
initialI=1; 
% Number of spread simulations to perform
numSims=1000; 
% Number of sampling simulations to perform
numRuns=1000; 
% Maximum run time for spread and sampling simulations
tFinal=5000;  
% Specify whether to display progress messages
progress = "no";

%% ------------------------------------------------------------------------
% END USER INPUT

%% ------------------------------------------------------------------------ 
% GENERATE SPREAD DATA

% Spread data for random sampling
simDataRandom = runSpread_SI_1(popSize,initialI,beta,numSims,tFinal,progress);
% Spread data for repeated sampling
simDataRepeated = runSpread_SI_1_REPEATED(popSize,initialI,beta,numSims,tFinal,progress);

%% ------------------------------------------------------------------------
% VARY SAMPLE SIZE

sampleInterval = 180; % 30, 90, 180 used in Supplementary Fig S15
sampleSizeVec = 10:10:100;
EDI_random = zeros(1,length(sampleSizeVec));
EDI_repeated = zeros(1,length(sampleSizeVec));

% Perfom sampling simulations
for i=1:length(sampleSizeVec)    
    % Random sampling
    [~,EDI,~,~,~] = runSampling_SI_1(simDataRandom,numRuns,popSize,sampleSizeVec(i),sampleInterval,tFinal,progress);
    EDI_random(i) = EDI;
    % Repeated sampling
    [~,EDI,~,~,~] = runSampling_SI_1_REPEATED(simDataRepeated,numRuns,popSize,sampleSizeVec(i),sampleInterval,tFinal,progress);
    EDI_repeated(i) = EDI;
end

% Plot
figure(); hold on; box off; grid on; set(gca,'Fontsize',16,'Linewidth',2);

% Define x plot values
xplot = sampleSizeVec;

% Plot EDP_random and EDP_repeated
myplot(1) = plot(xplot,100*EDI_random/popSize);
myplot(2) = plot(xplot,100*EDI_repeated/popSize);

% Define plot markers, colours and styles
myplot(1).Marker = 'o';
myplot(2).Marker = '^'; 

IBMblue = [100,143,255]/256;
myplot(1).Color = IBMblue;
myplot(1).MarkerFaceColor = IBMblue;
myplot(1).MarkerEdgeColor = IBMblue;

IBMpink = [221,37,128]/256;
myplot(2).Color = IBMpink;
myplot(2).MarkerFaceColor = IBMpink;
myplot(2).MarkerEdgeColor = IBMpink;

for l=1:2;  myplot(l).LineWidth = 2; myplot(l).LineStyle = ':'; myplot(l).MarkerSize = 10; end

% Title
title(['Sample interval =' 32 num2str(sampleInterval)]);

% Format x axis
xlabel('Sample size');
xlim([sampleSizeVec(1) sampleSizeVec(end)]);

% Format y axis
ylabel('EDP (% of total population)');

% Legend
legEntries = {'Random sampling','Repeated sampling'};
leg = legend([myplot(1) myplot(2)],legEntries);
leg.Location = 'northeast';

%% ------------------------------------------------------------------------
% VARY SAMPLE INTERVAL

sampleSize = 50; % 10, 20, 50 used in Supplementary Fig S15
sampleIntervalVec = 30:15:180;
EDI_random = zeros(1,length(sampleIntervalVec));
EDI_repeated = zeros(1,length(sampleIntervalVec));

% Perform sampling simulations
for i=1:length(sampleIntervalVec) 
    % Random sampling
    [~,EDI,~,~,~] = runSampling_SI_1(simDataRandom,numRuns,popSize,sampleSize,sampleIntervalVec(i),tFinal,progress);
    EDI_random(i) = EDI;
    % Repeated sampling
    [~,EDI,~,~,~] = runSampling_SI_1_REPEATED(simDataRepeated,numRuns,popSize,sampleSize,sampleIntervalVec(i),tFinal,progress);
    EDI_repeated(i) = EDI;
end

% Plot
figure(); hold on; box off; grid on; set(gca,'Fontsize',16,'Linewidth',2);

% Define x plot values
xplot = sampleIntervalVec;

% Plot EDP_random and EDP_repeated
myplot(1) = plot(xplot,100*EDI_random/popSize);
myplot(2) = plot(xplot,100*EDI_repeated/popSize);

% Define plot markers, colours and styles
myplot(1).Marker = 'o';
myplot(2).Marker = '^'; 

IBMblue = [100,143,255]/256;
myplot(1).Color = IBMblue;
myplot(1).MarkerFaceColor = IBMblue;
myplot(1).MarkerEdgeColor = IBMblue;

IBMpink = [221,37,128]/256;
myplot(2).Color = IBMpink;
myplot(2).MarkerFaceColor = IBMpink;
myplot(2).MarkerEdgeColor = IBMpink;

for l=1:2;  myplot(l).LineWidth = 2; myplot(l).LineStyle = ':'; myplot(l).MarkerSize = 10; end

% Title
title(['Sample size =' 32 num2str(sampleSize)]);

% Format x axis
xlabel('Sample interval');
xlim([sampleIntervalVec(1) sampleIntervalVec(end)]);
xticks = 30:15:180;
set(gca,'XTick', xticks);

% Format y axis
ylabel('EDP (% of total population)');

% Legend
legEntries = {'Random sampling','Repeated sampling'};
leg = legend([myplot(1) myplot(2)],legEntries);
leg.Location = 'southeast';
