function [c] = MPEC_constraint(x0,auxdata) 

%%%%%%%%%%%%%%%
% Constraints for the random coefficients Logit estimated via MPEC.
%%%%%%%%%%%%%%%


K=auxdata{3};
prods=auxdata{4};
T=auxdata{5};
numProdsTotal=auxdata{6};
x=auxdata{7};
IV=auxdata{8};
v=auxdata{9};
nn=auxdata{10};
share=auxdata{11};
prodsMarket=auxdata{12};
marketStarts=auxdata{13};
marketEnds=auxdata{14};
oo=auxdata{15};
sharesum=auxdata{16};
marketForProducts=auxdata{17};

theta1 = x0(1:K+1, 1);                              % mean tastes
theta2 = x0(K+2:2*K+2, 1);                          % st. deviation of tastes
delta = x0(2*K+3:2*K+2+numProdsTotal, 1);                 % mean utilities
g = x0(2*K+3+numProdsTotal:end, 1);                       % moment condition values

cong = g - IV'*(delta - x*theta1);  % constraints on moment conditions

expmu = exp(x*diag(theta2)*v);      % exponentiated deviations from mean utilities
expmeanval = exp(delta);
[EstShare, simShare] = ind_shnormMPEC(expmeanval,expmu, oo, sharesum, marketForProducts);


c = [EstShare - share ;
     cong ]; 
 
 
end