%% ------------------------------------------------------------------------
% DEFINE OBJECTIVE FUNCTION TO BE MINIMISED BY BAYESOPT

function OUTPUT = objFunc(simData,numRuns,popSize,numSentinels,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,~,baselineEDP)
    [~, ~, ~, ~, ECDP, ~] = runSampling_2(simData,numRuns,popSize,numSentinels,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,"no");
    OUTPUT = 100*(ECDP-baselineEDP)/baselineEDP;
end