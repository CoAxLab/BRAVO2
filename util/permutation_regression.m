function [p true_stat stat_ci perm] = permutation_regression(x,y,contrast,varargin);

% function [p true_stat CI boot] = permutation_regression(X,Y,C,OPTS);
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% PERMUTATION_REGRESSION:
% Premutes X against Y using OLS regression to evalaute the contrast
% effect modeled in the contrast vector C. 
% 
% INPUTS:
%        X    = NxC matrix where N = # observations and C = # of independent variables;
%        Y    = Nx1 vector.
%        C    = Cx1 contrast vector for statistical estimation
%
%     Optional Input: 
%           'scramble_cols' = which colums to permute from X
%           'niter'     = number of simulations to perform
%           'stat_type' = what kind of contrast to perform?  A 
%           t-test 't' or a simple beta comparison 'simple' (Default)
% 
% OUTPUTS:
%     p         = Significance of the observed relationship between X and Y
%
%     true_stat = Observed contrast statistics
%
%     CI        = Confidence interval of the bootstrap
%
%     perm      = niter x 1 vector of simulated contrasts
%     
% Note: Permutation happens by scrambling the order of each column
% of X independently.  Use the 'scramble_cols' vector to isolate specific
% columns to permute.
% 
% Written by Tim Verstynen (2007; updated 2013)
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.

niter = 1000;
scramble_cols = find(contrast);
stat_type = 'simple';

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

switch stat_type
  case 't';
    true_stat = contrast * [betas ./ se];
  case 'simple'
    true_stat = contrast * betas;
  otherwise
    error('Unknown comparison type');
end;

% Estimate 
true_se   = contrast * se;
sse  = sum(r.^2);
sst = sum(y.^2);
true_rsqr = 1-(sse/sst);


for iter = 1:niter
    
    nx = [];
    
    for xcol = 1:length(scramble_cols);
        nx = [nx x(randperm(size(x,1)),scramble_cols(xcol))];
    end;
    
    x(:,scramble_cols) = nx;  
    
    [betas, r] = ols_regress(y,x);
    se = se_calc(r,x)';

    sse  = sum(r.^2); 
    boot_rsqr(iter) = 1-(sse/sst);
    
    switch stat_type
      case 't';
        perm(iter) = contrast * [betas ./ se];

      case 'simple'
    	perm(iter) = contrast * betas;

        otherwise
        error('Unknown comparison type');
    end;

    boot_se(iter)   = contrast * se;
    
end;

if true_stat > mean(perm);
    p = length(find(perm > true_stat))/niter;
else
    p = length(find(perm < true_stat))/niter;
end;


% Estmiated CI to account for asymmetries in the simulations
sboot = sort(perm);
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

