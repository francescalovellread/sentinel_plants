%% ------------------------------------------------------------------------
% DEFINE DETERMINISTIC CONSTRAINT FUNCTION FOR BAYESOPT

function oC = objConstraint(X)
    oC = X.Ssamp<=X.Stot;
end
