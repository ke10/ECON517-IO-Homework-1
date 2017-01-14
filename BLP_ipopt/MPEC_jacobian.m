function J = MPEC_jacobian(x0, auxdata) 

%%%%%%%%%%%%%%
% Jacobian for the random coefficients Logit estimated via MPEC.
%%%%%%%%%%%%%%


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
 

    nx0 = size(x0,1);
    ng = size(g,1);
    ooo = ones(1,K+1);
    
    % Evaluate the Gradients
    dSdtheta2 = zeros(numProdsTotal,K+1);
    dSddeltaDIAG = zeros(numProdsTotal,prods);
    dSddelta = zeros(numProdsTotal, numProdsTotal);
    
    for t=1:T,
        index = marketStarts(t):marketEnds(t);
        ooo1 = ones(prodsMarket(t),1);
        for rr = 1:nn,
            dSddeltaDIAG(index,1:prodsMarket(t)) = dSddeltaDIAG(index,1:prodsMarket(t)) + (diag(simShare(index,rr)) - simShare(index,rr)*simShare(index,rr)')/nn;
            dSdtheta2(index,:) = dSdtheta2(index,:) + (simShare(index,rr)*ooo).*(ooo1*v(:,rr)').*( x(index,:) - (ooo1*(simShare(index,rr)'*x(index,:))))/nn;
        end
        dSddelta(index,index) = dSddeltaDIAG(index,1:prodsMarket(t));
    end


    dc1 = [zeros(numProdsTotal,K+1), dSdtheta2, dSddelta, zeros(numProdsTotal,ng)];
    dc2 = [IV'*x, zeros(ng, K+1), -IV', eye(ng)];
    J = [dc1; dc2];

    J=sparse(J);
    
    
end