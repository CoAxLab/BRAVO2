function [ci, p] = ci_percent(coefs, boot, alpha);

%function [CI, p_value] = ci_percentile(coefs, boot, alpha)
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% CI_PERCENTILE:
% Estimates the 95% confidence intervals and p-value using the discretely 
% sampled permutation vector.
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

% Set default value to 0.05 with two tails
if nargin < 3 | isempty(alpha); alpha = 0.025; end;

% First calculate the p-value
[B, n] = size(boot');
pct = sum(boot' < repmat(coefs,1,B)')./B;
p = pct;
is_lower = find(pct > 0.5);
p(is_lower) = 1-p(is_lower);

% Next the asymmetric CI
sboot = sort(boot')';
lb = round(size(sboot,2)*alpha);
if lb == 0; lb = 1; end;

ub = round(size(sboot,2)*(1-alpha));
ci = [sboot(:,lb)'; sboot(:,ub)'];

