%% Code for Fig 2: Schematic illustrating computational implementation of model
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 4th August 2022

%% -----------------------------------------------------------------------------------------------
% This code produces the epidemic curves and prevalence distribution shown in Fig 2.
clear; close all;

%% ------------------------------------------------------------------------
% DEFINE MODEL PARAMETERS

% Crop population size (PC)
Pc = 1000; 
% Sentinel population size (PS)
Ps = 100; 
% Compute total population size
popSize = Pc+Ps;
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

% Crop sample size
cropSampleSize = 50;
% Sentinel sample size
sentinelSampleSize = 50;
% Sample interval
sampleInterval = 100;

% Number of spread simulations to perform
numSims = 10000; 
% Number of sampling simulations to perform
numRuns = 10000;
% Maximum run time for spread and sampling simulations
tFinal = 5000; 
% Specify whether to display progress messages
progress = "no";

% Plotting parameters for incidence histogram
first = 0;
last = 50;
binWidth = 2;

%% ------------------------------------------------------------------------
% PERFORM SPREAD AND SAMPLING SIMULATIONS

simData = runSpread_SCI_2(popSize,Ps,D0,U0,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,numSims,tFinal,progress);
[sampleData,EDI,~,~,~,~,~] = runSampling_SCI_2(simData,numRuns,popSize,Ps,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,progress);


% Choose incidence curve to plot at random
selectedRuns = randsample(numSims,1);
selectedSimData = simData(selectedRuns,:);

% Extract time vector for plotting
xplot = selectedSimData{1,1};

%% ------------------------------------------------------------------------ 
% PLOT FIG 2A

figure();
hold on; box off; set(gca,'Fontsize',16,'Linewidth',2);

% Plot curves for Detectable and Undetectable crops and sentinels
Dc_plot = plot(xplot/365,selectedSimData{1,2});
Ds_plot = plot(xplot/365,selectedSimData{1,3});
Uc_plot = plot(xplot/365,selectedSimData{1,4});
Us_plot = plot(xplot/365,selectedSimData{1,5});

% Define plot colours and line widths
Dc_plot.Color = [100,143,255]/256; Dc_plot.LineWidth = 2;
Ds_plot.Color = [221,37,128]/256; Ds_plot.LineWidth = 2;
Uc_plot.Color = [120,93,241]/256; Uc_plot.LineWidth = 2;
Us_plot.Color = [254,97,0]/256; Us_plot.LineWidth = 2;

% Format x axis
xlabel('Time (years)');
xlim([0 10]);

% Format y axis
ylabel('Number of infected plants')
ylim([0 1200]);

% Legend
legEntries = {'','','',''};
leg = legend([Dc_plot,Uc_plot,Ds_plot,Us_plot],legEntries);
leg.Location = 'east';
leg.Box = "off";

%% ------------------------------------------------------------------------
% PLOT FIG 2B

figure();
hold on; box off; set(gca,'Fontsize',16,'Linewidth',2);

% Plot curves for total detectable plants and total infected crops
Det_plot = plot(xplot,selectedSimData{1,2}+selectedSimData{1,3});
Crop_plot = plot(xplot,selectedSimData{1,2}+selectedSimData{1,4});

% Define plot colours and line widths
Det_plot.Color = [0 0.5 0.9]; Det_plot.LineWidth = 2;
Crop_plot.Color = [0.5 0.5 0.5]; Crop_plot.LineWidth = 2; Crop_plot.LineStyle = "-.";

% Format x axis
xlabel('Time (days)');
xlim([0 600]);

% Format y axis
ylabel('Number of infected plants')
yticks = 0:25:150;
set(gca,'YTick', yticks);

% Legend
legEntries = {'',''};
leg = legend([Det_plot,Crop_plot],legEntries);
leg.Location = 'east';
leg.Box = "off";

%% ------------------------------------------------------------------------
% PLOT FIG 2C

figure();
hold on; box off; set(gca,'Fontsize',16,'Linewidth',2);

% Sort simulated detection prevalences ready to plot as a histogram
discoveryIncidences = 100*sampleData(:,1)/Pc; % Extract column vector of discovery incidences from sampleData
binEdges = first:binWidth:last+binWidth; % Define bin edges for histogram (last+binWidth ensures last is actually contained within the final bin; otherwise the bins may end before last.
incBinMidpoints = first+binWidth/2:binWidth:last+binWidth/2; % Compute bin midpoints
incidenceBars = histcounts(discoveryIncidences,binEdges); % Sort discovery incidences into bins
incidenceBars = incidenceBars/sum(incidenceBars); % Normalise

% Plot histogram of detection prevalences
% incidencePlot = bar(incBinMidpoints,incidenceBars);
incidencePlot = bar(mybar.XData,mybar.YData);
incidencePlot.FaceColor = [146 207 79]/256;
incidencePlot.EdgeColor = [0 0 0];
incidencePlot.LineWidth = 1;
incidencePlot.BarWidth = 1;

% Plot EDP
EDPplot = xline(100*EDI/Pc,'k-.','linewidth',4);
EDPplot.Alpha = 1;

% Format x axisl
xlabel('Detection prevalence (% of crop population)');
ylabel('Probability density');
xlim([first, last]);

% Legend
leg = legend([incidencePlot EDPplot],'Simulated detection prevalences','Expected detection prevalance (EDP)');
leg.Location = 'northeast'; leg.Box = 'off'; leg.FontSize = 16;
