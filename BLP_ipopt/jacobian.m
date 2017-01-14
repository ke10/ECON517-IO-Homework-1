function Ddelta = jacobian(theta2, x, v, delta, K, numProdsTotal, prods, T, nn, prodsMarket, marketStarts, marketEnds, oo, sharesum, marketForProducts)

%%%%%%%%%%%%%%
% Evaluates the Jacobian for the computation of standard errors
%%%%%%%%%%%%%%

expmu = exp(x*diag(theta2)*v);      % exponentiated deviations from mean utilities
expmeanval = exp(delta);
[EstShare, simShare] = ind_shnormMPEC(expmeanval,expmu, oo, sharesum, marketForProducts);

ooo = ones(1,K+1);

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


Ddelta = -inv(dSddelta)*dSdtheta2;



