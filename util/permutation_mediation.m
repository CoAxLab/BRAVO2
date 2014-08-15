function [coeffs, ab_perm, a_perm, b_perm, c_perm] = permutation_mediation(X,Y,M,C,niters)

%function [coeffs ab_boot, a_boot, b_boot, c_boot] = permutation_mediation(X,Y,M,C,niters)
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% PERMUTATION_MEDIATION:
% Performs a permutation mediation between two series of data
% to estimate the signficiance of mediator variables M1-Mi on the relationship
% between X & Y. Follows methods reported in Cerin et al. (2006) & Preacher
% and Hayes (2008). 
% 
% INPUTS:
%       X,Y,M  = independent, dependent and mediator variables
%       respectively.  X & Y are Nx1 vectors of continous data. M is an NxI
%       vector with I = # of mediator variables
%
%       C = NxC Matrix of covariates.  If no covariates desired,
%       then give an empty matrix (i.e., [])
%      
%       niters  = Number of permutations to run in the bootstrap
%       (Default 500 iterations).
%
% OUTPUT:  
%       
%       coeffs = Object of path coefficients a, ab, b, c, * cprime
%
%       ab_perm = permutation array of ab pathways
%
%       a_perm = permutation array of a pathways
%
%       b_perm = permutation array of b pathways
%
%       c_perm = permutation array of c-prime pathways
% 
% Written by T. Verstynen & A. Weinstein (2011). Updated 2013;
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.


if nargin < 5 | isempty(niters)
  niters = 1000;
end;

% Make column vectors
if size(X,2)>size(X,1); X = X'; end;
if size(M,2)>size(M,1); M = M'; end;
if size(Y,2)>size(Y,1); Y = Y'; end;
if size(C,2)>size(C,1); C = C'; end;

n_mediators = size(M,2);
n_covs      = size(C,2);

warning off
% estimate up X->M
[acoef, r] = ols_regress(M,[X C]);
a = acoef(1,:)'; Sa = se_calc(r,X)';

% estimate up X->Y
[ccoef, r] = ols_regress(Y,[X C]);
c = ccoef(1); Sc = se_calc(r,X);

% estimate up XM->Y
[prime_coef, r] = ols_regress(Y,[X M C]);
cprime = prime_coef(1); Scprime = se_calc(r, X);
b = prime_coef(2:n_mediators+1); Sb = se_calc(r,M); %prime_stats.se(2);

% store the paths
coeffs.a  = a(:);
coeffs.ab = a(:).*b(:);
coeffs.b  = b(:);
coeffs.c  = c;
coeffs.cprime = cprime;

% Setup the output variable
ab_perm = zeros(n_mediators,niters);
a_perm = zeros(n_mediators,niters);
b_perm = zeros(n_mediators,niters);
c_perm = zeros(1,niters);

% Now run the bootstrap
for it = 1:niters
  
    nx = X(randperm(length(X)));
    nm = M(randperm(size(M,1)),:);
    ny = Y(randperm(length(Y)));

    % X->M
    [acoef, r] = ols_regress(nm,[nx C]);
    a = acoef(1,:); 

    % XM->Y
    [prime_coef, r] = ols_regress(ny,[nx nm C]);
    c = prime_coef(1);
    b = prime_coef(2:n_mediators+1); 
    
    % Store each path separately
    ab_perm(:,it) = a(:).*b(:);
    a_perm(:,it) = b(:);
    b_perm(:,it) = a(:);
    c_perm(it)   = c;
    
end;
warning on

return;

function se = se_calc(resids, x);
% Calculates the standard error of each regression coefficient
n = length(x);
sigma = sqrt( nansum(resids.^2) / (n-2) );
se = sqrt( sigma.^2 ./ nansum(x.^2) );
return

function [betas, resid] = ols_regress(y,x);

ind = find(~isnan(sum([y x],2)));

% Calculates the OLS regression with multiple DVs
betas = inv(x(ind,:)'*x(ind,:))*x(ind,:)'*y(ind,:); 
y_hat = x(ind,:)*betas;

resid = NaN(size(y));
resid(ind,:) = y(ind,:)-y_hat;

return;


  
