%% Code for Figs 4BCD,5ABC: the optimal number of sentinels to include in the sample and the corresponding reductions in EDP
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Thursday 17th February 2022

%% -----------------------------------------------------------------------------------------------
% For a given number of sentinels added to the population, this code iterates over
% specified values of the sample size and sample interval and in each instance applies a
% Bayesian optimisation algorithm to compute the optimal number of sentinels to include in
% the sample and the corresponding reduction in EDP compared to the baseline level.

%% -----------------------------------------------------------------------------------------------
% READ IN DATA

% Provide path to data folder
cd '/Users/francescalovell-read/Documents/DPhil/Sentinels_paper/Code/new_results/Fig_4BCD_5ABC_results/Ps=0';

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
Pc = params.Pc;
Ps = params.Ps;

%% -----------------------------------------------------------------------------------------------
% EXTRACT DATA

% Extract sample size and sample interval vectors
sampleSizeVec = unique(dataStore(:,1))';
sampleIntervalVec = unique(dataStore(:,2))';

nSize = length(sampleSizeVec);
nInterval = length(sampleIntervalVec);

% Extract matrix of the optimal total number of sentinels in the sample
optTotSentinels = reshape(dataStore(:,3),[nInterval nSize]);

% Extract matrix of the reduction in EDP at the optimum
EDPreduction = reshape(dataStore(:,4),[nInterval nSize]);

% Extract matrix of the resultant EDP at the optimum
EDPresultant = reshape(dataStore(:,5),[nInterval nSize]);

%% ------------------------------------------------------------------------
% PLOT

% Define x plot values
xplot = sampleSizeVec;

% Define colours for IBM colourblind safe palette
IBMblue = [100,143,255]/256;
IBMpurple = [120,93,241]/256;
IBMpink = [221,37,128]/256;
IBMorange = [254,97,0]/256;
IBMyellow = [255,176,0]/256;

%% ------------------------------------------------------------------------
% PLOT OPTIMAL NUMBER OF SENTINELS TO INCLUDE IN SAMPLE

figure(); grid on; hold on;

% Grey out background in non-feasible region
if Ps==50
    pgon = polyshape([25 50 200 200 0],[25 50 50 200 200]); % For Ps=50
elseif Ps==100
    pgon = polyshape([25 100 200 200 0],[25 100 100 200 200]); % For Ps=100
elseif Ps==200
    pgon = polyshape([25 200 0],[25 200 200]); % For Ps=200
end

plotpg = plot(pgon);
plotpg.FaceColor = [0 0 0];
plotpg.FaceAlpha = 0.1;

% Plot line marking upper bound on the number of sentinels sampled
boundline1 = line([sampleSizeVec(1),Ps],[sampleSizeVec(1),Ps]);
boundline2 = line([Ps,sampleSizeVec(end)],[Ps,Ps]);
boundline1.LineWidth = 2.5; boundline2.LineWidth = 2.5;
boundline1.Color = [0 0 0]; boundline2.Color = [0 0 0];

% Plot optimal nuber of sentinels to include over range of sample sizes for each sample
% interval considered
myplot(1) = plot(xplot,optTotSentinels(1,:));
myplot(2) = plot(xplot,optTotSentinels(2,:));
myplot(3) = plot(xplot,optTotSentinels(3,:));

% Title
title(['P_S =' 32 num2str(Ps) 32 'sentinels']);

% Format x axis
xlabel('Sample size (N)');
xlim([sampleSizeVec(1) sampleSizeVec(end)]);
xticks = 25:25:200;
set(gca,'XTick', xticks, 'Fontsize', 16);

% Format y axis
ylabel('Optimal number of sentinels in sample (N_S^*)');
ylim([0 200]);
yticks = 0:25:200;
set(gca, 'YTick', yticks, 'Fontsize', 16);

% Legend
legEntries = {'min(P_S,N)','\Delta = 30 days','\Delta = 60 days','\Delta \geq 90 days'};
leg = legend([boundline1 myplot(1) myplot(2) myplot(3)],legEntries);
leg.Location = 'northwest';

% Define plot markers, colours and styles
myplot(1).Marker = '^';
myplot(2).Marker = '*'; 
myplot(3).Marker = 's';

myplot(1).Color = IBMblue;
myplot(1).MarkerFaceColor = IBMblue;
myplot(1).MarkerEdgeColor = IBMblue;

myplot(2).Color = IBMpink;
myplot(2).MarkerFaceColor = IBMpink;
myplot(2).MarkerEdgeColor = IBMpink;

myplot(3).Color = IBMyellow;
myplot(3).MarkerFaceColor = IBMyellow;
myplot(3).MarkerEdgeColor = IBMyellow;

for l=1:3;  myplot(l).LineWidth = 2; myplot(l).LineStyle = ':'; myplot(l).MarkerSize = 10; end

box off; set(gca,'Fontsize',16,'Linewidth',2);

%% ------------------------------------------------------------------------
% PLOT CORRESPONDING REDUCTIONS IN EDP

figure(); grid on; hold on;

% Plot baseline level
baseline = yline(0); baseline.LineWidth = 2; baseline.LineStyle = '--';

% Plot EDP reductions at optimum
myplot(1) = plot(xplot,EDPreduction(1,:));
myplot(2) = plot(xplot,EDPreduction(2,:));
myplot(3) = plot(xplot,EDPreduction(3,:));
myplot(4) = plot(xplot,EDPreduction(4,:));
myplot(5) = plot(xplot,EDPreduction(5,:));

% Title
title(['P_S =' 32 num2str(Ps) 32 'sentinels']);

% Format x axis
xlabel('Sample size (N)');
xlim([sampleSizeVec(1) sampleSizeVec(end)]);
xticks = 25:25:200;
set(gca,'XTick', xticks, 'Fontsize', 16);

% Format y axis
ylabel('Change in EDP from baseline at optimum (%)');

% Legend
legEntries = {'\Delta = 30 days','\Delta = 60 days','\Delta = 90 days','\Delta = 120 days','\Delta = 150 days'};
leg = legend([myplot(1) myplot(2) myplot(3) myplot(4) myplot(5)],legEntries);
leg.Location = 'southeast';

% Define plot markers, colours and styles
myplot(1).Marker = '^';
myplot(2).Marker = '*'; 
myplot(3).Marker = 's';
myplot(4).Marker = 'o'; 
myplot(5).Marker = 'x'; 

myplot(1).Color = IBMblue;
myplot(1).MarkerFaceColor = IBMblue;
myplot(1).MarkerEdgeColor = IBMblue;

myplot(2).Color = IBMpink;
myplot(2).MarkerFaceColor = IBMpink;
myplot(2).MarkerEdgeColor = IBMpink;

myplot(3).Color = IBMyellow;
myplot(3).MarkerFaceColor = IBMyellow;
myplot(3).MarkerEdgeColor = IBMyellow;

myplot(4).Color = IBMorange;
myplot(4).MarkerFaceColor = [1 1 1];
myplot(4).MarkerEdgeColor = IBMorange;

myplot(5).Color = IBMpurple;
myplot(5).MarkerFaceColor = IBMpurple;
myplot(5).MarkerEdgeColor = IBMpurple;

for l=1:5;  myplot(l).LineWidth = 2; myplot(l).LineStyle = ':'; myplot(l).MarkerSize = 10; end

box off; set(gca,'Fontsize',16,'Linewidth',2);

%% ------------------------------------------------------------------------
% PLOT RESULTANT EDP

figure(); grid on; hold on;

% Plot EDP reductions at optimum
myplot(1) = plot(xplot,100*EDPresultant(1,:)/Pc);
myplot(2) = plot(xplot,100*EDPresultant(2,:)/Pc);
myplot(3) = plot(xplot,100*EDPresultant(3,:)/Pc);
myplot(4) = plot(xplot,100*EDPresultant(4,:)/Pc);
myplot(5) = plot(xplot,100*EDPresultant(5,:)/Pc);

% Title
title(['P_S =' 32 num2str(Ps) 32 'sentinels']);

% Format x axis
xlabel('Sample size (N)');
xlim([sampleSizeVec(1) sampleSizeVec(end)]);
xticks = 25:25:200;
set(gca,'XTick', xticks, 'Fontsize', 16);

% Format y axis
ylabel('EDP at optimum (% of crop population)');

% Legend
legEntries = {'\Delta = 30 days','\Delta = 60 days','\Delta = 90 days','\Delta = 120 days','\Delta = 150 days'};
leg = legend([myplot(1) myplot(2) myplot(3) myplot(4) myplot(5)],legEntries);
leg.Location = 'northeast';

% Define plot markers, colours and styles
myplot(1).Marker = '^';
myplot(2).Marker = '*'; 
myplot(3).Marker = 's';
myplot(4).Marker = 'o'; 
myplot(5).Marker = 'x'; 

myplot(1).Color = IBMblue;
myplot(1).MarkerFaceColor = IBMblue;
myplot(1).MarkerEdgeColor = IBMblue;

myplot(2).Color = IBMpink;
myplot(2).MarkerFaceColor = IBMpink;
myplot(2).MarkerEdgeColor = IBMpink;

myplot(3).Color = IBMyellow;
myplot(3).MarkerFaceColor = IBMyellow;
myplot(3).MarkerEdgeColor = IBMyellow;

myplot(4).Color = IBMorange;
myplot(4).MarkerFaceColor = [1 1 1];
myplot(4).MarkerEdgeColor = IBMorange;

myplot(5).Color = IBMpurple;
myplot(5).MarkerFaceColor = IBMpurple;
myplot(5).MarkerEdgeColor = IBMpurple;

for l=1:5;  myplot(l).LineWidth = 2; myplot(l).LineStyle = ':'; myplot(l).MarkerSize = 10; end

box off; set(gca,'Fontsize',16,'Linewidth',2);

%% ------------------------------------------------------------------------
% WRITE RESULTS TO TXT FILES

cd '/Users/francescalovell-read/Documents/DPhil/Sentinels_paper/Code/new_results/Fig_4BCD_5ABC_results/all';

writematrix(sampleSizeVec,'sampleSizeVec.txt','delimiter','\t');
writematrix(sampleIntervalVec,'sampleIntervalVec.txt','delimiter','\t');
writematrix(EDPreduction,['valueAtOptStore',num2str(Ps),'.txt']);
