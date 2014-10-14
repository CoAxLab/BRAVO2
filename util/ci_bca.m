function [ci, p] = ci_bca(coefs, boot, alpha);

%function [CI, p_value] = ci_bca(coefs, boot, alpha)
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% CI_BCA:
% Estimates the confidence intervals and p-value according to the
% bias corrected and accelerated formula described in DiCiccio & Efron
% (1996).
%
% INPUTS:
%       coefs  = observed path values [Nx1]
%
%       boot   = bootstrap distributions [NxI]
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
% Written by T. Verstynen & A. Weinstein (2011).
%
% All code is licensed under the GNU Public License (version 3.0).  

n_coef = length(coefs);

if ~coefs & n_coef == 1 & size(boot,2)>1
    coefs = repmat(coefs,1,size(boot,2));
end;

% Set default value to 0.05 with two tails
if nargin < 3 | isempty(alpha); alpha = 0.025; end;

% Modified from Tor Wager routine

% Bias factor
[B, ncols] = size(boot);
prop_less  = sum(boot < repmat(mean(boot),B,1))./B;
z0 = norminv(prop_less);

% Acceleration Factor: Borrowed from bootbca_pval.m by Tor Wager et al.
jstat = jackknife(@mean, boot);
n = size(jstat,1);
score = -(jstat - repmat(mean(jstat), n, 1) ); % score function at stat;
skew = sum(score.^3)./(sum(score.^2).^1.5);  %skewness of the score function
a =  skew ./ 6;  % acceleration

% Next estimate the adjusted CI
z(1) = norminv(alpha);
z(2) = norminv(1 - alpha);

a1 = z0 + ( (z0 + z(1)) ./ (1 - a .* (z0 + z(1))) );
a1 = normcdf(a1);

a2 = z0 + ( (z0 + z(2)) ./ (1 - a .* (z0 + z(2))) );
a2 = normcdf(a2);

ci = [diag(prctile(boot,a1*100))'; diag(prctile(boot,a2*100))']; 

% Then estimate the p-values
pct = sum(boot < repmat(coefs,n,1))./B;

% Adjust for ceilings
pct_lb = max(pct,1./B);
pct_ub = min(pct,1 - 1./B);

% Correct for ceiling/floor hits
if ~pct_lb; pct_lb = 0.0001; end;
if ~pct_ub; pct_ub = 0.0001; end;
if pct_lb==1; pct_lb = 0.9999; end;
if pct_ub==1; pct_ub = 0.9999; end;

% Now adjust
zpct = norminv(pct_ub) - z0;
zadj = ( zpct .* (1- a.* (z0)) - z0 ) ./ (1 + a.* zpct);

p = normcdf(zadj);
is_lower = find(pct > 0.5);
p(is_lower) = 1-p(is_lower);

