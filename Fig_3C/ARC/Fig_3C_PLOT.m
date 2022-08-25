%% Code for Fig 3C: computing the baseline EDP for a range of sample sizes and intervals
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 4th August 2022

%% -----------------------------------------------------------------------------------------------
% This code takes in the output from 'Fig_3C_RUN.m' and uses it to plot Fig 3C

%% -----------------------------------------------------------------------------------------------
% READ IN DATA

% Provide path to data folder
cd '/Users/francescalovell-read/Documents/DPhil/Sentinels_paper/Code/new_results/Fig_3C_results/';

% Make directory of all results files
FileList  = dir(fullfile(pwd,'results_*.txt'));

% Import and store data from results files
dataStore = [];
for ID = 1:numel(FileList)
    filename = sprintf('results_%d.txt',ID);
    data = importdata(filename);
    dataStore = [dataStore;data]; 
end

% Read in parameter table
params = readtable('params.txt');
P = params.P;

%% -----------------------------------------------------------------------------------------------
% EXTRACT DATA

% Extract sample size and sample interval vectors
sampleSizeVec = unique(dataStore(:,1))';
sampleIntervalVec = unique(dataStore(:,2))';

nSize = length(sampleSizeVec);
nInterval = length(sampleIntervalVec);

% Extract and smooth EDP matrix
EDPmatrix = reshape(dataStore(:,3),[nInterval nSize]);
EDPmatrix_smooth = imgaussfilt(EDPmatrix,1);

%% ------------------------------------------------------------------------
% PLOT

figure();

% Create contour plot of EDPs expressed as percentage of total population
levels = [0 1 2:2:30];
[C,h] = contourf(100*EDPmatrix_smooth/P,levels);

% Set contour and label styles and specify manual labelling
h.LineColor = [0 0 0]; h.LineStyle = '-'; h.LineWidth = 1;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end

% Define colourbar
colbar = colorbar;
ylabel(colbar, 'Baseline EDP (% of population)', 'Fontsize', 16);
set(colbar,'linewidth',2,'fontsize',16);

% Format x axis
xlabel('Sample size (N)');
xticks = 1:5:nSize;
xticklabels = strsplit(num2str(sampleSizeVec(xticks)));
set(gca,'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

% Format y axis
ylabel('Sample interval (\Delta days)');
yticks = 1:3:nInterval;
yticklabels = strsplit(num2str(sampleIntervalVec(yticks)));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'Fontsize', 16);

box off; set(gca,'Fontsize',16,'Linewidth',2);
