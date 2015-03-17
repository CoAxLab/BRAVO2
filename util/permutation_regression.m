function [p true_stat stat_ci perm] = permutation_regression(x,y,contrast,varargin);

% Emulate changes made in bootstrap_regession

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
%               (Note: Do not include an intercept column in X)
%        Y    = Nx1 vector.
%        C    = Cx1 contrast vector for statistical estimation
%
%     Optional Input: 
%           'n_iter'     = number of simulations to perform
%           'stat_type' = what kind of contrast to perform?  A 
%                         t-test 't' or a simple beta comparison 
%                         'simple' (Default)
%           'reg_type'  = type of regression to use: 'ols_regress' (simple OLS, Default)
%                          or 'qr_regress' (QR decomposition)
% 
% OUTPUTS:
%     p         = Significance of the observed relationship between X and Y
%
%     true_stat = Observed contrast statistics
%
%     CI        = Confidence interval of the bootstrap
%
%     perm      = n_iter x 1 vector of simulated contrasts
%     
% Note: Permutation happens by scrambling the order of each column
% of X independently.  Use the 'scramble_cols' vector to isolate specific
% columns to permute.
% 
% Written by Tim Verstynen (2007; updated 2013)
%
% Revised and released as BRAVO 2.0 by T. Verstynen (2014)
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.

% Reset random number genrator as a precaution (changing to using Jason's approach)
RandStream('mt19937ar','Seed',sum(100*clock));

n_iter = 1000;
stat_type = 'simple';
reg_type = 'ols_regress';
n_threads = 0; % Note multithreading is not fully supported yet. Don't turn this on.
is_par = 0;

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

switch stat_type
  case 't';
    true_stat = [0 contrast] * [betas ./ se];
  case 'simple'
    true_stat = [0 contrast] * betas;
  otherwise
    error('Unknown comparison type');
end;

% Estimate 
true_se   = contrast * se;

try 
  if ~matlabpool('size') & n_threads > 0;
    if n_threads == Inf;
      matlabpool OPEN;
    else
      matlabpool(n_threads);
    end;
  end
  is_par = 1;
catch
  is_par = 0;
end;

qr_flag = strcmpi(reg_type,'qr_regress');

if is_par 
  parfor iter = 1:n_iter
     [perm(iter), perm_se(iter)] = simulate_iteration(x,y,contrast,qr_flag);
  end;

else
  for iter = 1:n_iter
    [perm(iter), perm_se(iter)] = simulate_iteration(x,y,contrast,qr_flag);
  end;
end;

switch stat_type
  case 't';
    perm = perm./perm_se;
end;

if true_stat > 0;
    p = length(find(perm < true_stat))/n_iter;
else
    p = length(find(perm > true_stat))/n_iter;
end;

% Estmiated CI to account for asymmetries in the simulations
sboot = sort(perm);
lb_pt = round(length(sboot)*0.025);
ub_pt = round(length(sboot)*0.975);

stat_ci = [sboot(lb_pt) sboot(ub_pt)];
    
return;

% -----------------------------------------
function [boot, boot_se] = simulate_iteration(x,y,contrast,qr)

  nx = [];
  for xcol = 1:size(x,2)
    nx = [nx x(randperm(size(x,1)),xcol)];
  end;

  ny = y(randperm(size(y,1)));

  switch qr
  case 1
    [betas,r] = qr_regress(ny,nx);
  case 0
    [betas,r] = ols_regress(ny,nx);
  end;

  boot = [0 contrast] * betas;

  % Calculate SE
  n = size(nx,1);
  sigma = sqrt( nansum(r.^2) / (n-2) );

  for i = 1:size(nx,2);
    se(i) = sqrt( sigma.^2 ./ nansum(nx(:,i).^2) );
  end;

  boot_se   = contrast * se';


