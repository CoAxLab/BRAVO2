function [betas, resid] = ols_regress(y,x);

% function [betas, resid] = ols_regress(y,x);
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% OLS_REGRESS:
% Ordinary least squares sub-function. This is a core utility function 
% for the mediation, regression and correlation routines.
% 
% Released as BRAVO 2.0 by T. Verstynen (2014)
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.

% Insert intercept
x = x2fx(x);

% Filter for empty data
ind = find(~isnan(sum([y x],2)));

% Calculates the OLS regression with multiple DVs
betas = inv(x(ind,:)'*x(ind,:))*x(ind,:)'*y(ind,:); 
y_hat = x(ind,:)*betas;
resid = NaN(size(y));
resid(ind,:) = y(ind,:)-y_hat;