%%%%%%%%%%%%%%%%%%%%%%%
% MAIN SCRIPT
%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

%%%%%%%%%%% 
% SETTINGS
%%%%%%%%%%%

randn('seed',155)                                   % reset normal random number generator
rand('seed',155)                                    % reset uniform random number generator
nn = 100;                                           % # draws to simulate shares
tol_inner = 1.e-14;                                 % NFP inner loop tolerance, for starting values here                                            % # Markets
prods = 25;  
T = 50;
starts = 5;                                         % # random start values to check during estimation


%%%%%%%%%%%%%%%%%%%%%
% Number of products per market
% Discrete uniform between 1 and prods
%%%%%%%%%%%%%%%%%%%%%%

prodsMarket = randi([1 prods],T,1);  % Number of products per market,


% Use two dimensional matrices for speed reasons (vectorization)
% Market indices tell us the first and last spots of records for each
% market in list of all products
marketStarts = zeros(T,1);
marketEnds = zeros(T,1);
marketStarts(1) = 1;
marketEnds(1) = prodsMarket(1);
for t=2:T
    marketStarts(t) = marketEnds(t-1) + 1;
    marketEnds(t) = marketStarts(t) + prodsMarket(t) - 1;
end;
numProdsTotal = marketEnds(T);  % Total number of products in all markets

%%%%%%%%%%%%%%%%%%%%%%%%
% TRUE PARAMETER VALUES
%%%%%%%%%%%%%%%%%%%%%%%%
costparams = ones(6,1)/2;
betatrue = [2 1.5 1.5 .5 -3]';                     % true mean tastes
betatrue0 = betatrue;
K = size(betatrue,1)-1;                             % # attributes
covrc = diag( [ .5 .5 .5 .5 .2] );                  % true co-variances in tastes
rctrue = covrc(find(covrc));
thetatrue = [betatrue;sqrt(rctrue)];
v = randn(length(betatrue),nn);                     % draws for share integrals during estimation
rc = chol(covrc)'*v;                                % draws for share integrals for data creation
sigmaxi = 0.5;
covX = -.8*ones(K-1)+1.8*eye(K-1);                  % allow for correlation in X's
covX(1,3) = .3;
covX(2,3) = .3;
covX(3,1) = .3;
covX(3,2) = .3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create indices to speed share calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oo = ones(1,nn);                                    % index for use in simulation of shares

sharesum = sparse(zeros(T,numProdsTotal));  % used to create denominators in logit predicted shares (i.e. sums numerators)
for t=1:T
    sharesum(t,marketStarts(t):marketEnds(t)) = 1;
end

marketForProducts = zeros(numProdsTotal,1);  % market indices for each product
% used to expand one logit denominator per consumer to one for each product
for t=1:T
    marketForProducts(marketStarts(t):marketEnds(t)) = t;
end


%%%%%%%%%%%%%%%%%%%%%
% SIMULATE DATA
%%%%%%%%%%%%%%%%%%%%%
randn('seed',5000)                                   % reset normal random number generator
rand('seed',5000)                                    % reset uniform random number generator
xi = randn( numProdsTotal,1)*sigmaxi;                % draw demand shocks
A = randn(numProdsTotal,K-1)*chol(covX);             % product attributes
prand = rand(numProdsTotal,1)*5;
price = 3 +   xi*1.5 +  prand + sum(A,2);
z = rand(numProdsTotal,length(costparams)) + 1/4*repmat( abs( prand +  sum(A,2)*1.1 ) ,1,length(costparams));
x = [ones(numProdsTotal,1) A price];
[share,nopurch] = mksharesim(betatrue, x, xi, rc, prods, T, sharesum, marketForProducts);
y = log(share) - log(nopurch(marketForProducts,:) );     % log-odds ratio for no RC logit estimation
% matrix of instruments
iv = [ones(numProdsTotal,1) A z A.^2 A.^3 z.^2 z.^3 prod(A,2) prod(z,2) kron(A(:,1),ones(1,size(z,2))).*z  kron(A(:,2),ones(1,size(z,2))).*z];
IV = iv;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2SLS ESTIMATION OF HOMOGENEOUS LOGIT MODEL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pz = iv*inv(iv'*iv)*iv';
xhat = pz*x;
PX2sls = inv(xhat'*xhat)*xhat';                     % project. matrix on weighted x-space for 2SLS
beta2sls = PX2sls*y;
se2sls = sqrt(diag( mean((y-x*beta2sls).^2)*inv(x'*pz*x) ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choices for GMM estimation
% GMM weighting matrix follows Nevo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expmeanval0 = exp(y);     % starting value for mean values 
W = inv(IV'*IV);  % GMM weighting matrix
PX = inv(x'*IV*W*IV'*x)*x'*IV*W*IV'; % Used to calculate starting values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPEC ESTIMATION OF RANDOM COEFFICIENTS LOGIT MODEL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta20 = 0.5*abs(beta2sls);
startvalues = repmat(theta20',starts,1).* cat(2,ones(size(theta20)),rand([starts-1,size(theta20',2)] )' *1 )';
nIV = size(IV,2);       % # instrumental variables

GMPEC = 1.0e20;
CPUtMPEC = 0;
FuncEvalMPEC = 0;
GradEvalMPEC = 0;
HessEvalMPEC = 0;
expmeanval = expmeanval0;
X_MPEC = zeros(2*(K+1)+ numProdsTotal + nIV, starts);
% logshare = log(share);

for reps = 1:starts
    %%% NOTE: START VALUES SET TO GIVE SAME INITIAL GMM OBJECTIVE FUNCTIONAL
    %%% VALUE AS WITH NFP
    theta20 = startvalues(reps,:)';  % starting values for standard deviations
    delta0 = invertshares(theta20, x, expmeanval, tol_inner, share, v, oo, sharesum, marketForProducts);  % starting values for mean utilities using BLP inversion
    theta10 = PX*delta0; % starting values for linear parameters using IV like regression of delta on covariates
    resid0 = delta0 - x*theta10; % starting values for xi terms, not directly used
    g0 = IV'*resid0;  % starting values for GMM moments
    x0 = [theta10; theta20; delta0; g0];  % master starting values
    
    x_L = -Inf*ones(length(x0),1);   % Lower bounds for x.
    x_L(K+2:2*K+2) = 0;              % standard deviations are nonnegative
    x_U =  Inf*ones(length(x0),1);   % Upper bounds for x.
    
    
    
    nx0 = size(x0,1);   % number of MPEC optimization parameters
    
    % Sparsity patterns of constraints, 0 if derivative 0, 1 otherwise
    % Derivatives of market shares
    
    c11 = zeros(numProdsTotal, K+1);   % market shares with respect to mean parameters
    c12 = ones(numProdsTotal, K+1);  % market shares with respect to standard deviations
    
    c13 = zeros(numProdsTotal,numProdsTotal); % market shares with respect to mean utilities
    for t=1:T
        c13(marketStarts(t):marketEnds(t),marketStarts(t):marketEnds(t)) = 1;
    end
    
    c14 = zeros(numProdsTotal, nIV); % market shares with respect to moment values
    
    % Derivatives of moments
    
    c21 = ones(nIV, K+1);   % moments with respect to mean parameters 
    c22 = zeros(nIV, K+1);    % moments with respect to standard deviations
    c23 = ones(nIV, numProdsTotal);  % moments with respect to mean utilities
    c24 = eye(nIV);       % moments with respect to moment values 
    
    ConsPattern = [ c11 c12 c13 c14; c21 c22 c23 c24 ];   
        
    % Hessian pattern
    
    HessianPattern = zeros(nx0, nx0);
    HessianPattern(2*K+3+numProdsTotal:end, 2*K+3+numProdsTotal:end) = ones(nIV, nIV);
    HessianPattern(K+2:2*K+2,K+2:2*K+2+numProdsTotal) = ones(K+1, K+1+numProdsTotal);
    HessianPattern(2*K+3:2*K+2+numProdsTotal,K+2:2*K+2) = ones(numProdsTotal, K+1);
    for t=1:T,
        range = marketStarts(t):marketEnds(t);
        index = 2*K+2+range;
        HessianPattern(index, index) = ones(prodsMarket(t), prodsMarket(t));
    end

    
%%%%%%%%%%%%%%%%    
% OPTIMIZATION %
%%%%%%%%%%%%%%%%

% Bounds on x
options.lb = x_L;    % Lower Bounds for x
options.ub = x_U;    % Upper Bounds for x
  
% Bounds on constraints
options.cl = zeros(size(ConsPattern,1),1);             % Lower bounds on constraints.
options.cu = zeros(size(ConsPattern,1),1);             % Upper bounds on constraints.

% Set IPOPT options.
options.ipopt.jac_c_constant        = 'no';
options.ipopt.hessian_constant      = 'no';
options.ipopt.mu_strategy           = 'adaptive';
options.ipopt.derivative_test       = 'none';         % 'second-order' 'none'
options.ipopt.tol                   = 1e-6;
options.ipopt.max_iter              = 100; 
auxdata = { ConsPattern W K prods T numProdsTotal x IV v nn share prodsMarket marketStarts marketEnds oo sharesum marketForProducts HessianPattern};
options.auxdata = auxdata;
  
% Callback functions
funcs.objective = @MPEC_objective;                       % Objective Function
funcs.constraints = @MPEC_constraint;                    % Constraints
funcs.gradient = @MPEC_gradient;                         % Gradient
funcs.jacobian = @MPEC_jacobian;                         % Jacobian
funcs.jacobianstructure = @MPEC_jacobian_structure;      % Structure of Jacobian
funcs.hessian = @MPEC_hessian;                           % Hessian
funcs.hessianstructure  = @MPEC_hessian_structure;       % Structure of Hessian
  
% Run Optimization
[X, info] = ipopt(x0,funcs,options);
    
X_MPEC(:, reps) = X;
objective(reps,1) = MPEC_objective(X, auxdata);
status(reps,1) = info.status;

theta1MPEC_rep = X(1:K+1);
theta2MPEC_rep = X(K+2:2*K+2);
deltaMPEC_rep = X(2*K+3:2*K+2+numProdsTotal, 1); 
GMPEC_rep = MPEC_objective(X, auxdata);
INFOMPEC_rep = info.status;

if (GMPEC_rep < GMPEC && INFOMPEC_rep==0),
     thetaMPEC1 = theta1MPEC_rep;
     thetaMPEC2 = theta2MPEC_rep;
     deltaMPEC = deltaMPEC_rep;
     GMPEC = GMPEC_rep;
     INFOMPEC = INFOMPEC_rep;
end

end     % End Loop over different startvalues


% COVARIANCE MATRIX FOR MPEC STRUCTURAL PARAMETERS
deltaSE = invertshares(thetaMPEC2, x, expmeanval, tol_inner, share, v, oo, sharesum, marketForProducts);
resid = deltaSE - x*thetaMPEC1;  % xi
delta = deltaSE;
Ddelta = jacobian(thetaMPEC2, x, v, delta, K, numProdsTotal, prods, T, nn, prodsMarket, marketStarts, marketEnds, oo, sharesum, marketForProducts);
covg = zeros(size(IV,2));
for ii =1:length(IV),
    covg = covg + IV(ii,:)'*IV(ii,:)*(resid(ii)^2);
end
Dg = [x Ddelta]'*IV;            % gradients of moment conditions
covMPEC = inv( Dg*W*Dg')*Dg*W*covg*W*Dg'*inv( Dg*W*Dg');



  
