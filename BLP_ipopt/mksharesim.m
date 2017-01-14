function [share,nopurch] = mksharesim(betatrue,x,xi,rc, prods, T, sharesum, marketForProducts)

%%%%%%%%%%%%
% Generate market share data.
%%%%%%%%%%%%


delta = x*betatrue+xi;

MU = x*rc;
numer = exp( repmat(delta,1,size(rc,2)) + MU );  % exp of utility for each simulated consumer


sum1 = sharesum*numer;                  % sum of utility for each consumer
sum11 = 1./(1+sum1);                    % this is the denominator of the shares
denom1 = sum11(marketForProducts,:);          % this expands the denominator
shareForDraws = numer.*denom1;               % simulated shares for each draw
share = mean(shareForDraws,2);            % expected share (i.e. mean across draws)
nopurch = mean(sum11,2);           % Fraction in each market not purchasing


