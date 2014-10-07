function [p true_stat stat_ci boot] = bagging_regression(x,y,contrast,varargin);

% function [p true_stat CI sim] = bagging_regression(X,Y,C,OPTS);
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% BAGGING_REGRESSION:
% Uses an MCMC subsampling approach X against Y using OLS regression to evalaute the contrast
% effect modeled in the contrast vector C. 
% 
% INPUTS:
%        X    = NxC matrix where N = # observations and C = # of independent variables;
%               (Note: Do not include an intercept column in X)
%        Y    = Nx1 vector.
%        C    = Cx1 contrast vector for statistical estimation
%
%     Optional Input: 
%           'niter'     = number of simulations to perform
%           'stat_type' = what kind of contrast to perform?  A 
%                         t-test 't' or a simple beta comparison 
%                         'simple' (Default)
%           'ratio'     = subsampling ratio (Default = 2/3rds)
%           'reg_type'  = type of regression to use: 'ols_regress' (simple OLS, Default)
%                          or 'qr_regress' (QR decomposition)
% 
% OUTPUTS:
%     p         = Significance of the observed relationship between X and Y
%
%     true_stat = mean of bootstrapped statistic;
%
%     CI        = Confidence interval of the bootstrap
%
%     sim      = niter x 1 vector of simulated contrasts
%     
% Note: Bootstrapping happen by subsampling with replacement the input vectors
% to estimate a confidence interval of a parameter.
% 
% Written by Tim Verstynen (2013)
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.

% Reset random number genrator as a precaution (changing to using Jason's approach)
RandStream('mt19937ar','Seed',sum(100*clock));

niter = 1000;
stat_type = 'simple';
ratio = 2/3;
reg_type = 'ols_regress';

% Get variable input parameters
for v=1:2:length(varargin),
    eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
end

if ~sum(strcmpi(reg_type,{'ols_regress','qr_regress'}));
    error('Unknown regression type. Options are ols_regress and qr_regress');
end;

% Check to make sure you're dealing with column data;
if size(x,2) > size(x,1); x = x'; end;
if size(y,2) > size(y,1); y = y'; end;

% Get the initial model
eval(sprintf('[betas, r] = %s(y,x);',reg_type));
se = se_calc(r,x)';

% Get number of data points
n = length(x);

for iter = 1:niter
    
    indx = randperm(n);
    indx = indx(1:round(ratio*n));

    nx = x(indx,:);
    ny = y(indx,:);

    eval(sprintf('[betas, r] = %s(ny,nx);',reg_type));
    se = se_calc(r,nx)';

    switch stat_type
      case 't';
        boot(iter) = [0 contrast] * [betas ./ se];

      case 'simple'
    	boot(iter) = [0 contrast] * betas;

        otherwise
        error('Unknown comparison type');
    end;

    boot_se(iter)   = contrast * se;
end;

if mean(boot) > 0;
    p = length(find(boot < 0))/niter;
else
    p = length(find(boot > 0))/niter;
end;

true_stat = mean(boot);

% Estmiated CI to account for asymmetries in the simulations
sboot = sort(boot);
lb_pt = round(length(sboot)*0.025);
ub_pt = round(length(sboot)*0.975);

stat_ci = [sboot(lb_pt) sboot(ub_pt)];

    
return;




