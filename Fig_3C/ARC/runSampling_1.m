function [sampleData, EDP, EDT, EDR, perc95] = runSampling_1(simData,numRuns,popSize,sampleSize,sampleInterval,tFinal,progress)

% INPUT
% simData: cell array of simulated incidence curves, generated from runSpread_1.m
% numRuns: number of sampling runs to perform (must be <= number of simulations available
% in simData!)
% popSize: total population size
% sampleSize: number of plants to sample on each sampling round
% sampleInterval: sampling interval (time between each sampling round)
% tFinal: final permissible sample time
% progress: specifies whether progress messages are displayed ("yes" or "no")

% OUTPUT
% sampleData: a three-column matrix in which each row corresponds to a single
% sampling run. Entries in the first column are (total) discovery prevalences; entries in the
% second column are the corresponding discovery times; entries in the third column are the
% numbers of sampling rounds that occurred.
% EDP: the simulated expected discovery prevalence
% EDT: the simulated expected discovery time
% EDR: the simulated expected number of sampling rounds

rng('shuffle');

if (progress ~= "yes" && progress ~= "no")
    fprintf('ERROR: Please enter a valid argument for progress ("yes" or "no")\n\n'); return
end

timerSample1 = tic;

N = sampleSize; Delta = sampleInterval; P = popSize;
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
% uniformly distributed on the interval [0,Delta] to mimic pathogen
% introduction at a random time relative to the sampling scheme.
initSampleTimes = rand(numRuns,1)*Delta;

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
        numDetectable = selectedSimData{i,2}(sampleIndex); % Compute total number of 'Detectable' plants
        numUndetectable = selectedSimData{i,3}(sampleIndex); % Compute total number of 'Undetectable' plants
        stateVec = ones(1,numDetectable); stateVec(P) = 0; % Create vector containing a 1 for every 'Detectable' plant and a 0 for every 'Undetectable'/'Healthy' plant
        selectVec = randsample(stateVec,N); % Take a random sample of size N (without replacement) from this vector
        detSample = sum(selectVec); % Total number of 'Detectable' plants in the sample
        if detSample > 0 % Then a 'Detectable' plant has been sampled
            detectionIndicator = 1;
        else
            sampleTime = sampleTime + Delta; % Move on to next sample time
            samplingRound = samplingRound + 1;
        end
    end
    sampleData(i,1) = numDetectable + numUndetectable;
    sampleData(i,2) = sampleTime;
    sampleData(i,3) = samplingRound;
end

EDP = mean(sampleData(:,1)); % Expected total discovery prevalence ('Detectable' and 'Undetectable')
EDT = mean(sampleData(:,2)); % Expected discovery time
EDR = mean(sampleData(:,3)); % Expected number of sampling rounds
perc95 = prctile(sampleData(:,1),95);

elapsedTimeSample1 = toc(timerSample1);
if progress == "yes"
    fprintf(strcat('DONE! (',num2str(elapsedTimeSample1),32,'secs)\n\n'));
end
end