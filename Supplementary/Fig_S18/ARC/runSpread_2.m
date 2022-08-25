function simData = runSpread_2(popSize,numSentinels,initialD,initialU,betaC,betaS,epsilonC,epsilonS,gammaC,gammaS,numSims,tFinal,progress)
rng('shuffle');
% INPUT
% popSize: total population size.
% numSentinels: number of sentinels in the population
% initialD: total initial number of 'Detectable' individuals (crops and sentinels)
% initialU: total inital number of 'Undetectable' individuals (crops and sentinels)
% betaC: transmission coefficient for 'Detectable' crops
% betaS: transmission coefficient for 'Detectable' sentinels
% epsilonC: transmission coefficient scaling factor for 'Undetectable' crops
% epsilonS: transmission coefficient scaling factor for 'Undetectable' sentinels
% gammaC: 'Undetectable' period for crops
% gammaS: 'Undetectable' period for sentinels
% numSims: number of simulations to run
% tFinal: maximum time for each simulation to run
% progress: specifies whether progress messages are displayed ("yes" or "no")

% OUTPUT
% sim_data: a (no_sims x 7) cell array containing simulation results. Each row corresponds
% to a run; columns contain vectors for t, Dc ('Detectable' crops), Ds ('Detectable' sentinels),
% Uc ('Undetectable' crops), Us ('Undetectable' sentinels), Hc ('Healthy' crops) and Hs
% ('Healthy' sentinels).

rng('shuffle');

if (progress ~= "yes" && progress ~= "no")
    fprintf('ERROR: Please enter a valid argument for progress ("yes" or "no")\n\n'); return
end

timerSpread2 = tic;

P = popSize; Ps = numSentinels; Pc = P-Ps; sentinel_prop = Ps/P;
bc = betaC; bs = betaS; ec = epsilonC; es = epsilonS; gc = gammaC; gs = gammaS;
% Create cell array for storing results
simData = cell(numSims,7);
if progress == "yes"
    fprintf('Running spread simulations...\t')
end

% Run stochatic simulations
for i=1:numSims
    % On each run, divide initial infected individuals between crops and sentinels
    % according to their relative proportions in the population
    D0 = initialD; U0 = initialU;
    Dc0 = 0; Ds0 = 0; Uc0 = 0; Us0 = 0; randsD = rand(1,D0); randsU = rand(1,U0);
    for j = 1:D0
        if randsD(j)<=sentinel_prop
            Ds0 = Ds0+1;
        else
            Dc0 = Dc0+1;
        end
    end
    for j = 1:U0
        if randsU(j)<=sentinel_prop
            Us0 = Us0+1;
        else
            Uc0 = Uc0+1;
        end
    end
    % Compute initial numbers of susceptible crops and sentinels
    Hc0 = Pc-Dc0-Uc0; Hs0 = Ps-Ds0-Us0;
    
    t = 0; Dc = Dc0; Ds = Ds0; Uc = Uc0; Us = Us0; Hc = Hc0; Hs = Hs0;
    vecLength = 1+2*(Hc0+Hs0)+Uc0+Us0;
    tvec = zeros(1,vecLength);
    Dcvec = zeros(1,vecLength); Dsvec = zeros(1,vecLength);
    Ucvec = zeros(1,vecLength); Usvec = zeros(1,vecLength);
    Hcvec = zeros(1,vecLength); Hsvec = zeros(1,vecLength);
    index = 1; tvec(index)=t; Dcvec(index)=Dc0; Dsvec(index)=Ds0; Ucvec(index)=Uc0; Usvec(index)=Us0; Hcvec(index)=Hc0; Hsvec(index)=Hs0;
   
    while (t<tFinal && Hc+Hs+Uc+Us>0)
        % Compute individual reaction propensities
        a1 = Hc*(bc*Dc+bs*Ds+ec*bc*Uc+es*bs*Us);
        a2 = Hs*(bc*Dc+bs*Ds+ec*bc*Uc+es*bs*Us);
        a3 = (1/gc)*Uc;
        a4 = (1/gs)*Us;
        a0 = a1+a2+a3+a4; % Compute total reaction propensity
        r1 = rand(1); tau = (1/a0)*log(1/r1); % Compute time to next reaction
        r2 = rand(1);
        if r2<=a1/a0 % Then a susceptible crop becomes a cryptic crop
            Hc = Hc-1; Uc = Uc+1;
        elseif r2<=(a1+a2)/a0 % Then a susceptible sentinel becomes a cryptic sentinel
            Hs = Hs-1; Us = Us+1;
        elseif r2<=(a1+a2+a3)/a0 % Then a cryptic crop becomes a symptomatic crop
            Uc = Uc-1; Dc = Dc+1;
        else % A cryptic sentinel becomes a symptomatic sentinel
            Us = Us-1; Ds = Ds+1;
        end
        t = t+tau;
        index = index + 1;
        tvec(index) = t;
        Dcvec(index) = Dc; Dsvec(index) = Ds;
        Hcvec(index) = Hc; Hsvec(index) = Hs;
        Ucvec(index) = Uc; Usvec(index) = Us;
    end
    simData{i,1} = tvec;
    simData{i,2} = Dcvec; simData{i,3} = Dsvec;
    simData{i,4} = Ucvec; simData{i,5} = Usvec;
    simData{i,6} = Hcvec; simData{i,7} = Hsvec;
end

elapsedTimeSpread2 = toc(timerSpread2);
if progress == "yes"
    fprintf(strcat('DONE! (',num2str(elapsedTimeSpread2),32,'secs)\n',num2str(numSims),32,'incidence curves generated.\n\n'));
end
end
