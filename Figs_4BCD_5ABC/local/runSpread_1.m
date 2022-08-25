function simData = runSpread_1(popSize,initialD,initialU,beta,epsilon,gamma,numSims,tFinal,progress)

% INPUT
% popSize: total population size
% initialD: initial number of 'Detectable' individuals
% initialU: initial number of 'Undetectable' individuals
% beta: transmission coefficient for 'Detectable' individuals
% epsilon: transmission coefficient scaling factor for 'Undetectable' individuals
% gamma: duration of 'Undetectable' period
% numSims: number of simulations to run
% tFinal: maximum time for each simulation to run
% progress: specifies whether progress messages are displayed ("yes" or "no")

% OUTPUT
% sim_data: a (no_sims x 4) cell array containing simulation results. Each
% row corresponds to a run; first, second, third and fourth columns contain
% vectors for t, D, U and H respectively.

rng('shuffle');

if (progress ~= "yes" && progress ~= "no")
    fprintf('ERROR: Please enter a valid argument for progress ("yes" or "no")\n\n'); return
end

timerSpread1 = tic;

P = popSize; D0 = initialD; U0 = initialU; b = beta; e = epsilon; g = gamma;
% Compute initial number of 'Healthy' individuals
H0 = P-D0-U0;
% Create cell array for storing results
simData = cell(numSims,4);
if progress == "yes"
    fprintf('Running spread simulations...\t')
end
% Run stochatic simulations
for i=1:numSims
    t = 0; D = D0; U = U0; H = H0;
    tvec = zeros(1,1+2*H0+U0); Dvec = zeros(1,1+2*H0+U0); Uvec = zeros(1,1+2*H0+U0); Hvec = zeros(1,1+2*H0+U0);
    index = 1; tvec(index)=t; Dvec(index)=D0; Uvec(index)=U0; Hvec(index)=H0;
    while (t<tFinal && (H+U)>0)
        a1 = b*H*D; a2 = b*e*H*U; a3 = (1/g)*U; % Compute individual reaction propensities
        a0 = a1+a2+a3; % Compute total reaction propensity
        r1 = rand(1); tau = (1/a0)*log(1/r1); % Compute time to next reaction
        r2 = rand(1);
        if r2<(a1+a2)/a0 % Then an H becomes a U
            U = U+1; H = H-1;
        else % Then a U becomes an D
            D = D+1; U = U-1;
        end
        t = t+tau;
        index = index + 1; 
        tvec(index) = t; 
        Dvec(index) = D; Uvec(index) = U; Hvec(index) = H; 
    end
    simData{i,1} = tvec;
    simData{i,2} = Dvec;
    simData{i,3} = Uvec;
    simData{i,4} = Hvec;
end
elapsedTimeSpread1 = toc(timerSpread1);
if progress == "yes"
    fprintf(strcat('DONE! (',num2str(elapsedTimeSpread1),32,'secs)\n',num2str(numSims),32,'incidence curves generated.\n\n'));
end
end
