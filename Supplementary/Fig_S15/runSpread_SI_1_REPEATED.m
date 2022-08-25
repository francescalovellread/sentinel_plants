function simData = runSpread_SI_1_REPEATED(popSize,initialI,beta,numSims,tFinal,progress)

% INPUT
% popSize: total population size.
% initialI: initial number of infected individuals
% beta: transmission coefficient
% numSims: number of simulations to run
% tFinal: maximum time for each simulation to run
% progress: specifies whether progress messages are displayed ("yes" or "no")

% OUTPUT
% sim_data: a (no_sims x 4) cell array containing simulation results. Each
% row corresponds to a run; first, second, third columns contain vectors
% for t, I and S respectively. Fourth column contains data on which individuals
% become infected at each time point.

if (progress ~= "yes" && progress ~= "no")
    fprintf('ERROR: Please enter a valid argument for progress ("yes" or "no")\n\n'); return
end

tic

P = popSize; I0 = initialI; b = beta;
% Compute initial number of susceptible individuals
S0 = P-I0;
% Create cell array for storing results
simData = cell(numSims,4);
if progress == "yes"
    fprintf('Running spread simulations...\t')
end
% Run stochatic simulations
for i=1:numSims
    t = 0; I = I0; S = S0;
    tvec = zeros(1,P+1-I0); Ivec = zeros(1,P+1-I0); Svec = zeros(1,P+1-I0); 
    storevec = cell(P+1-I0,1);
    % Select which plants begin as infected
    zerovec = zeros(1,P);
    choose = randsample(P,I0);
    zerovec(choose) = 1;
    storevec{1}=zerovec;
    
    index = 1; tvec(index)=t; Ivec(index)=I0; Svec(index)=S0;
    while (t<tFinal && S>0)
        a0 = b*S*I;
        r1 = rand(1);
        tau = (1/a0)*log(1/r1); 
        t = t+tau; I = I+1; S = S-1; 

        currentvec = storevec{index};
        uninf = find(~currentvec); 
        if length(uninf)==1
            choose = uninf; 
        else 
            choose = randsample(uninf,1); 
        end
        currentvec(choose)=1;
        newvec = currentvec;

        index = index + 1; tvec(index) = t; Ivec(index) = I; Svec(index) = S;
        storevec{index}=newvec;

    end
    simData{i,1} = tvec;
    simData{i,2} = Ivec;
    simData{i,3} = Svec;
    simData{i,4}=storevec;
end
elapsedTime = toc;
if progress == "yes"
    fprintf(strcat('DONE! (',num2str(elapsedTime),32,'secs)\n',num2str(numSims),32,'incidence curves generated.\n\n'));
end
end
