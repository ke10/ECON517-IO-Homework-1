function [EstShare, simShare] = ind_shnormMPEC(expmeanval,expmu, oo, sharesum, marketForProducts)

%%%%%%%%%%%%%%%
% Compute the individual probabilities of choosing each product.
% The probabilities are those associated with the normally-distributed r.c. logit.
%%%%%%%%%%%%%%%


numer = (expmeanval*oo ).*expmu;        % this is the numerator (oo speeds-up expanding mean utility by number of draws)
sum1 = sharesum*numer;
sum11 = 1./(1+sum1);                    % this is the denominator of the shares
denom1 = sum11(marketForProducts,:);          % this expands the denominator
simShare = numer.*denom1;               % simulated shares for each draw
EstShare = mean(simShare,2);            % expected share (i.e. mean across draws)
