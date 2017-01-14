function grad_f = MPEC_gradient(x0, auxdata)

%%%%%%%%%%%%%%
% Gradient for the random coefficients Logit esitmated via MPEC.
%%%%%%%%%%%%%%


W=auxdata{2};
K=auxdata{3};
numProdsTotal=auxdata{6};

theta1 = x0(1:K+1, 1);                      % mean tastes
theta2 = x0(K+2:2*K+2, 1);                  % st. deviation of tastes
delta = x0(2*K+3:2*K+2+numProdsTotal, 1);         % mean utilities
g = x0(2*K+3+numProdsTotal:end, 1);               % moment condition values

f = g'*W*g; 

    nx0 = size(x0,1);
    grad_f = zeros(nx0,1);
    grad_f(2*K+3+numProdsTotal:end,1) = 2*W*g;

end