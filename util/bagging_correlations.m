function [p true_cxy boot] = bagging_correlations(x,y,varargin);

%function [p r sim] = bagging_correlations(X,Y,OPTS);
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% BAGGING_CORRELATIONS:
% Performs a MCMC subsampling mediation between two series of data
% to estimate the signficiance of mediator variable M on the relationship
% between X & Y. Follows methods reported in Cerin et al. (2006) & Preacher
% and Hayes (2008). 
% 
% INPUTS:
%       x,y = Nx1 vectors to be correlated (both are permuted without
%       replacement)
%
%    Optional Input:
%       'niter' = number of permutations to run. Default = 1000;
%       'flag'   = 'spearmans','pearsons' (default)
%       'ratio'  = subsampling ratio (default = 67% of dataset)
%
% OUTPUT:
%       p  =  probability that the bootstrapped distribution crosses zero
%
%       r = true, unpermuted correlation value
%
%       sim = niter x 1 vector of simulated correlations
%
% Written by Timothy Verstynen (2013)
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.

% Reset random number genrator as a precaution (changing to using Jason's approach)
RandStream('mt19937ar','Seed',sum(100*clock));

% Defaults
ratio = 2/3;
flag = 'pearsons';
niter = 1000;

% Get variable input parameters
for v=1:2:length(varargin),
    eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
end

% Get rid of NaNs
good_vals = find(~isnan(x) & ~isnan(y));
x = x(good_vals);
y = y(good_vals);

n = length(x);

% Run the MCMC
boot =  [];
for it = 1:niter
  
  indx = randperm(n);
  indx = indx(1:round(ratio*n));

  nx = x(indx);
  ny = y(indx);
  
  switch flag
   case 'spearmans'
    cxy = spear(nx(:),ny(:));
    
   case 'pearsons'
    r = corrcoef(nx,ny);
    cxy = r(1,2);
  end;
  
  boot = [boot cxy];
end;

true_cxy = mean(boot);

% Calculate the probability of observing a stronger correlation by chance
if mean(boot) > 0;
    p = length(find(boot < 0))/niter;
else
    p = length(find(boot > 0))/niter;
end;
