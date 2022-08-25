%% ------------------------------------------------------------------------
% DEFINE OBJECTIVE FUNCTION TO BE MINIMISED BY BAYESOPT

function OUTPUT = objFunc(popSize,numSentinels,initialI,initialC,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,numSims,tFinal,progress,numRuns,cropSampleSize,sentinelSampleSize,sampleInterval,baselineEDP)
    if numSentinels==0
        OUTPUT = 0;
    else
        simData = runSpread_2(popSize,numSentinels,initialI,initialC,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,numSims,tFinal,progress);
        [~, ~, ~, ~, ECDP, ~] = runSampling_2(simData,numRuns,popSize,numSentinels,cropSampleSize,sentinelSampleSize,sampleInterval,tFinal,progress); 
        OUTPUT = 100*(ECDP-baselineEDP)/baselineEDP;
    end
end