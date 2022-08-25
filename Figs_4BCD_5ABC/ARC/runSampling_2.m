function [sampleData, EDP, EDT, EDR, ECDP, ESDP, perc95] = runSampling_2(simData,numRuns,popSize,numSentinels,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,progress)

% INPUT
% simData: cell array of simulated incidence curves, generated from runSpread_2.m
% numRuns: number of sampling runs to perform (must be <= number of simulations available in simData!)
% popSize: total population size
% numSentinels: number of sentinels in the population
% cropSampleSize: number of crop plants to sample on each sampling round
% sentinelSampleSize: number of sentinel plants to sample on each sampling round
% sampleInterval: sampling interval (time between each sampling round)
% tFinal: final permissible sample time
% progress: specifies whether progress messages are displayed ("yes" or "no")

% OUTPUT
% sampleData: a five-column matrix in which each row corresponds to a single sampling run.
% Entries in the first column are total discovery prevalences; entries in the second column are
% the corresponding discovery times; entries in the third column are the numbers of sampling
% rounds that occurred. Entries in columns four and five are the crop and sentinel
% discovery prevalences (i.e. total incidence broken down into species). All discovery
% prevalences include both 'Detectable' and 'Undetectable' individuals.
% EDP: the simulated expected total discovery prevalence (crops and sentinels, 'Detectable' and 'Undetectable')
% EDT: the simulated expected discovery time
% EDR: the simulated expected number of sampling rounds
% ECDP: the simulated expected crop discovery prevalence ('Detectable' and 'Undetectable')
% ESDP: the simulated expected sentinel discovery prevalence ('Detectable' and 'Undetectable'))

rng('shuffle');

if (progress ~= "yes" && progress ~= "no")
    fprintf('ERROR: Please enter a valid argument for progress ("yes" or "no")\n\n'); return
end

timerSample2 = tic;

Nc = cropSampleSize; Ns = sentinelSampleSize; Delta = sampleInterval;
P = popSize; Ps = numSentinels; Pc = P-Ps;
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
sampleData = zeros(numRuns,5);
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
       
        % Compute numbers of 'Detectable' and 'Undetectable' crops and sentinels at that time
        numDetectableCrops = selectedSimData{i,2}(sampleIndex); % Compute total number of 'Detectable' crop plants
        numDetectableSentinels = selectedSimData{i,3}(sampleIndex); % Compute total number of 'Detectable' sentinel plants
        numUndetectableCrops = selectedSimData{i,4}(sampleIndex); % Compute total number of 'Undetectable' crop plants
        numUndetectableSentinels = selectedSimData{i,5}(sampleIndex); % Compute total number of 'Undetectable' sentinel plants
       
        cropStateVec = ones(1,numDetectableCrops); cropStateVec(Pc) = 0; % Create vector containing a 1 for every 'Detectable' crop plant and a 0 for every 'Undetectable' crop plant
        sentinelStateVec = ones(1,numDetectableSentinels); sentinelStateVec(Ps) = 0; % Create vector containing a 1 for every 'Detectable' sentinel plant and a 0 for every 'Undetectable' sentinel plant
        
        % Take sample from crop plants
        selectCropVec = cropStateVec(randsample(length(cropStateVec),Nc)); % Take a random sample of size Nc (without replacement) from the crop vector
        detCropSample = sum(selectCropVec); % Total number of 'Detectable' crop plants in the sample
        
        % Take sample from setinel plants
        selectSentinelVec = sentinelStateVec(randsample(length(sentinelStateVec),Ns)); % Take a random sample of size Ns (without replacement) from the sentinel vector
        detSentinelSample = sum(selectSentinelVec); % Total number of 'Detectable' sentinel plants in the sample
             
        if detCropSample+detSentinelSample > 0 % Then an symptomatic plant has been sampled
            detectionIndicator = 1;
        else
            sampleTime = sampleTime + Delta; % Move on to next sample time
            samplingRound = samplingRound + 1;
        end
    end
    
    sampleData(i,1) = numDetectableCrops+numDetectableSentinels+numUndetectableCrops+numUndetectableSentinels; % Total number of infected individuals
    sampleData(i,2) = sampleTime; % Discovery time
    sampleData(i,3) = samplingRound; % Number of sampling rounds
    sampleData(i,4) = numDetectableCrops+numUndetectableCrops; % Total number of infected crops
    sampleData(i,5) = numDetectableSentinels+numUndetectableSentinels; % Total number of infected sentinels
end
EDP = mean(sampleData(:,1)); % Expected total discovery prevalence (crops and sentinels, 'Detectable' and 'Undetectable')
EDT = mean(sampleData(:,2)); % Expected discovery time
EDR = mean(sampleData(:,3)); % Expected number of sampling rounds
ECDP = mean(sampleData(:,4)); % Expected crop discovery prevalence ('Detectable' and 'Undetectable')
ESDP = mean(sampleData(:,5)); % Expected sentinel discovery prevalence ('Detectable' and 'Undetectable')
perc95 = prctile(sampleData(:,1),95);

elapsedTimeSample2 = toc(timerSample2);
if progress == "yes"
    fprintf(strcat('DONE! (',num2str(elapsedTimeSample2),32,'secs)\n',num2str(numRuns),32,'sampling simulations performed.\n\n'));
end
end