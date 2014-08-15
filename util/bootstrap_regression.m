function [p true_stat stat_ci boot] = bootstrap_regression(x,y,contrast,varargin);

% function [p true_stat CI boot] = bootsrap_regression(X,Y,C,OPTS);
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% BOOTSTRAP_REGRESSION:
% Bootstraps X against Y using OLS regression to evalaute the contrast
% effect modeled in the contrast vector C. 
% 
% INPUTS:
%        X    = NxC matrix where N = # observations and C = # of independent variables;
%        Y    = Nx1 vector.
%        C    = Cx1 contrast vector for statistical estimation
%
%     Optional Input: 
%           'niter'     = number of simulations to perform
%           'stat_type' = what kind of contrast to perform?  A 
%           t-test 't' or a simple beta comparison 'simple' (Default)
%           'ratio'     = subsampling ratio (Default = 2/3rds)
% 
% OUTPUTS:
%     p         = Significance of the observed relationship between X and Y
%
%     true_stat = mean of bootstrapped statistic;
%
%     CI        = Confidence interval of the bootstrap
%
%     boot      = niter x 1 vector of simulated contrasts
%     
% Note: Bootstrapping happen by subsampling with replacement the input vectors
% to estimate a confidence interval of a parameter.
% 
% Written by Tim Verstynen (2013)
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.

niter = 1000;
stat_type = 'simple';
ratio = 2/3;

% Check for appropriate input;
if nargin < 2
    error('Not enough inputs given');
end;

% Get variable input parameters
for v=1:2:length(varargin),
    eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
end

% Check to make sure you're dealing with column data;
if size(x,2) > size(x,1); x = x'; end;
if size(y,2) > size(y,1); y = y'; end;

% Get the initial model(OLS)
[betas, r] = ols_regress(y,x);
se = se_calc(r,x)';

% Get number of data points
n = length(x);

for iter = 1:niter
    
    indx = randperm(n);
    indx = indx(1:round(ratio*n));

    nx = x(indx,:);
    ny = y(indx,:);
    
    %[betas, dev, stats] = glmfit(x,y,'normal','constant','off');
    [betas, r] = ols_regress(ny,nx);
    se = se_calc(r,nx)';

    switch stat_type
      case 't';
        boot(iter) = contrast * [betas ./ se];

      case 'simple'
    	boot(iter) = contrast * betas;

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

%stat_ci = [mean(boot_stat)-std(boot_stat)*1.96 mean(boot_stat)+std(boot_stat)*1.96];
    
return;


function se = se_calc(resids, x);

if size(x,1)<size(x,2); x = x'; end;

n = size(x,1);
sigma = sqrt( nansum(resids.^2) / (n-2) );

for i = 1:size(x,2);
    se(i) = sqrt( sigma.^2 ./ nansum(x(:,i).^2) );
end;

return

function [betas, resid] = ols_regress(y,x);

ind = find(~isnan(sum([y x],2)));

% Calculates the OLS regression with multiple DVs
betas = inv(x(ind,:)'*x(ind,:))*x(ind,:)'*y(ind,:); 
y_hat = x(ind,:)*betas;

resid = NaN(size(y));
resid(ind,:) = y(ind,:)-y_hat;

return;

