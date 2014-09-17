function [ci, p] = ci_percentile(coefs, boot, alpha);

%function [CI, p_value] = ci_percentile(coefs, boot, alpha)
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% CI_PERCENTILE:
% Estimates the 95% confidence intervals and p-value according to the
% asymmetric CI estimates of Preacher & Hayes (2008) using a continuous 
% gaussian distribution. 
%
% INPUTS:
%       coefs  = observed path values [Nx1]
%
%       boot   = permutation distributions [NxI]
%
%       alpha  = alpha for CI estimation (default = 0.05);
%
% OUTPUT:
%
%       CI      = (1-alpha)% asymmetric confidence interval (lower, upper)
%
%       p_value = probability of simulating a stronger effect in the boot
%                 distribution than observed in the unpermuted data.
% 
% Written by T. Verstynen & A. Weinstein (2011). Updated 2013
%
% All code is licensed under the GNU Public License (version 3.0).  

% Set default value to 0.05 with two tails
if nargin < 3 | isempty(alpha); alpha = 0.025; end;

% First estimate the distribution parameters assuming 
mu = mean(boot);
sigma = std(boot);

p = normcdf(coefs,mu,sigma);
ci(1,:) = norminv(alpha,mu,sigma);
ci(2,:) = norminv(1-alpha,mu,sigma);

if coefs < 0; p = 1-p; end;
