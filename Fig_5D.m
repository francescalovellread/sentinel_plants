%% Code for Fig 5D: which of 50, 100 and 200 sentinels performs best
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 4th August 2022

%% -----------------------------------------------------------------------------------------------
% This code takes in the output from 'Figs_4BCD_5ABC.m' and uses it to plot Fig 5D
clear; close all;

%% -----------------------------------------------------------------------------------------------
% READ IN DATA

% Provide path to data folder
cd '/Users/francescalovell-read/Documents/DPhil/Sentinels_paper/Code/new_results/Fig_3BCD_4ABC_results/all';

% Read in sample size vector
sizeID = fopen('sampleSizeVec.txt','r');
sizeformatspec = '%d';
sampleSizeVec = fscanf(sizeID,sizeformatspec);
fclose(sizeID);
nSize = length(sampleSizeVec);

% Read in sample interval vector
intID = fopen('sampleIntervalVec.txt','r');
intformatspec = '%d';
sampleIntervalVec = fscanf(intID,intformatspec);
fclose(intID); 
nInterval = length(sampleIntervalVec);

% Read in valueAtOptStore matrices
valueAtOptStore50 = readmatrix('valueAtOptStore50.txt');
valueAtOptStore100 = readmatrix('valueAtOptStore100.txt');
valueAtOptStore200 = readmatrix('valueAtOptStore200.txt');

%% -----------------------------------------------------------------------------------------------
% DETERMINE WHICH NUMBER OF SENTINELS PERFORMS BEST AT EACH POINT

% Create empty matrix to store indicator values
sentinelIndicator = zeros(nInterval,nSize);
 
for i=1:nInterval % Iterate over sample intervals
    for j=1:nSize % Iterate over sample sizes
        
        EDP50 = valueAtOptStore50(i,j);
        EDP100 = valueAtOptStore100(i,j);
        EDP200 = valueAtOptStore200(i,j);
        
        % Determine which number of sentinels performed better and fill in storage matrix accordingly
        if min([EDP50 EDP100 EDP200]) == EDP50
            ind=1;
        elseif min([EDP50 EDP100 EDP200]) == EDP100
            ind=2;
        else
            ind=3;
        end
        sentinelIndicator(i,j)=ind;
    end
end

%% ------------------------------------------------------------------------
% PLOT

figure(); 

% Define colours for Tol colourblind safe palette
TolTeal = [96,168,153]/256;
TolBlue = [151,203,234]/256;
TolGreen = [15,119,51]/256;

% Define colormap
colormap([TolGreen;TolTeal;TolBlue]);

% Plot heat map of optimal number of sentinels 
clims=[1 3];
imagesc(flip(sentinelIndicator),clims); hold on

% Format x axis
xlabel('Sample size (N)')
xticks = 1:1:nSize;
xticklabels = strsplit(num2str(sampleSizeVec'));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

% Format y axis
ylabel('Sample interval (\Delta)');
yticks = 1:1:nInterval;
yticklabels = strsplit(num2str(sampleIntervalVec'));
set(gca, 'YTick', yticks, 'YTickLabel', flip(yticklabels), 'Fontsize', 16);

% Add text labels to each region
text1 = text(1.5,1,{'50 sentinels','gives greatest','reduction'});
text2 = text(3.5,2,{'100 sentinels gives','greatest reduction'});
text3 = text(6.5,3.5,{'200 sentinels gives','greatest reduction'});
text4 = text(7,5,{'50 sentinels','gives greatest','reduction'});

% Format text labels
text1.FontSize = 16; text1.HorizontalAlignment = 'center'; text1.Color = [1 1 1];
text2.FontSize = 16; text2.HorizontalAlignment = 'center';
text3.FontSize = 16; text3.HorizontalAlignment = 'center';
text4.FontSize = 16; text4.HorizontalAlignment = 'center'; text4.Color = [1 1 1];

box off; set(gca,'Fontsize',16,'Linewidth',2);