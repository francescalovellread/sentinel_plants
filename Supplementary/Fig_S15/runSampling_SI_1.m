function [sampleData, EDI, EDT, EDR, perc95] = runSampling_SI_1(simData,numRuns,popSize,sampleSize,sampleInterval,tFinal,progress)

% INPUT
% simData: cell array of simulated incidence curves, generated from runSpread_SI_1.m
% numRuns: number of sampling runs to perform (must be <= number of simulations available
% in simData!)
% popSize: total population size
% sampleSize: number of plants to sample on each sampling round
% sampleInterval: sampling interval (time between each sampling round)
% tFinal: final permissible sample time
% progress: specifies whether progress messages are displayed ("yes" or
% "no")

% OUTPUT
% sampleData: a three-column matrix in which each row corresponds to a single
% sampling run. Entries in the first column are discovery incidences; entries in the
% second column are the corresponding discovery times; entries in the third column are the
% numbers of sampling rounds that occurred.
% EDI: the simulated expected discovery incidence
% EDT: the simulated expected discovery time
% EDR: the simulated expected number of sampling rounds

if (progress ~= "yes" && progress ~= "no")
    fprintf('ERROR: Please enter a valid argument for progress ("yes" or "no")\n\n'); return
end

tic

N = sampleSize; D = sampleInterval; P = popSize;
% Extract number of simulations from simData
dataSize = size(simData); numSims = dataSize(1);
% If numRuns>numSims, return error message and exit
if numRuns > numSims
    fprintf(strcat('Error: number of sampling runs exceeds available number of simulations. Please set numRuns<=',num2str(numSims),'.\n\n'));
    return; 
end
selectedRuns = randsample(numSims,numRuns);
selectedSimData = simData(selectedRuns,:);
% Create empty matrix to store sampling data
sampleData = zeros(numRuns,3);
% Generate initial sampling times for each run. Initial sampling times are
% uniformly distributed on the interval [0,D] to mimic pathogen
% introduction at a random time relative to the sampling scheme.
initSampleTimes = rand(numRuns,1)*D;

if progress == "yes"
    fprintf('Running sampling simulations...\t')
end
% Run sampling simulations
for i=1:numRuns
    sampleTime = initSampleTimes(i);
    detectionIndicator = 0;
    samplingRound = 1;
    while (detectionIndicator == 0 && sampleTime <= tFinal)
        sampleIndex = sum(selectedSimData{i,1} <= sampleTime); % Determine which time point on the inidence curve corresponds to the sample time
        numInfected = selectedSimData{i,2}(sampleIndex); % Compute total number of infected plants
        stateVec = ones(1,numInfected); stateVec(P) = 0; % Create vector containing a 1 for every infected plant and a 0 for every healthy plant
        selectVec = randsample(stateVec,N); % Take a random sample of size N (without replacement) from this vector
        infSample = sum(selectVec); % Total number of infected plants in the sample
        if infSample > 0 % Then an infected plant has been sampled
            detectionIndicator = 1;
        else
            sampleTime = sampleTime + D; % Move on to next sample time
            samplingRound = samplingRound + 1;
        end
    end
    sampleData(i,1) = numInfected;
    sampleData(i,2) = sampleTime;
    sampleData(i,3) = samplingRound;
end
EDI = mean(sampleData(:,1)); % Expected discovery incidence
EDT = mean(sampleData(:,2)); % Expected discovery time
EDR = mean(sampleData(:,3)); % Expected number of sampling rounds
perc95 = prctile(sampleData(:,1),95);

elapsedTime = toc;
if progress == "yes"
    fprintf(strcat('DONE! (',num2str(elapsedTime),32,'secs)\n\n'));
end
end