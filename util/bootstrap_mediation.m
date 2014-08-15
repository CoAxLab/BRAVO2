function [coeffs, ab_boot, a_boot, b_boot, c_boot] = bootstrap_mediation(X,Y,M,C,niters,ratio)

%function [coeffs ab_boot, a_boot, b_boot, c_boot] = bootstrap_mediation(X,Y,M,C,niters,ratio)
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% BOOTSTRAP_MEDIATION:
% Performs a bootstrap mediation between two series of data
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
%        'ratio'     = subsampling ratio (Default = 2/3rds)
%
% OUTPUT:  
%       
%       coeffs = Object of path coefficients a, ab, b, * cprime
%
%       ab_boot = bootstrapped array of ab pathways
%
%       a_boot = bootstrapped array of a pathways
%
%       b_boot = bootstrapped array of b pathways\
%
%       c_boot = bootstrapped array of c-prime pathways
% 
% Written by T. Verstynen & A. Weinstein (2011). Modified by T. Verstynen (2013)
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.


%
% All code is licensed under the GNU Public License (version 3.0).  

% If no subsmaple ratio is given, go with 2/3rds
if nargin < 6
  ratio = 2/3;
end;

% If number of iteration is empty, return that
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

% Setup the output variable
ab_boot = zeros(n_mediators,niters);
a_boot = zeros(n_mediators,niters);
b_boot = zeros(n_mediators,niters);
c_boot = zeros(1,niters);

% Get number of samples;
n = size(X,1);

% Now run the bootstrap
for it = 1:niters

    %Subsampled index
    indx = randperm(n);
    indx = indx(1:round(ratio*n));

    % Pull from data
    nx = X(indx,:);
    nm = M(indx,:);
    ny = Y(indx,:);
    
    if ~isempty(C);
        nc = C(indx,:);
    else
        nc = [];
    end;
    
    % X->M
    [acoef, r] = ols_regress(nm,[nx nc]);
    a = acoef(1,:); 

    % XM->Y
    [prime_coef, r] = ols_regress(ny,[nx nm nc]);
    c = prime_coef(1);
    b = prime_coef(2:n_mediators+1); 
    
    % Store each path separately
    ab_boot(:,it) = a(:).*b(:);
    a_boot(:,it) = a(:);
    b_boot(:,it) = b(:);
    c_boot(it)   = c;

end;

% determine the mean coefficients
coeffs.a  = mean(a_boot,2);
coeffs.ab = mean(ab_boot,2);
coeffs.b  = mean(b_boot,2);
coeffs.cprime = mean(c_boot,2);

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


  
