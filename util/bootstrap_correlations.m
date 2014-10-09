function [p true_cxy boot] = bootstrap_correlations(x,y,varargin);

%function [p r sim] = bootstrap_correlations(X,Y,OPTS);
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% BOOTSTRAP_CORRELATIONS:
% Performs a bootstrapped mediation between two series of data
% to estimate the signficiance of mediator variable M on the relationship
% between X & Y. Follows methods reported in Cerin et al. (2006) & Preacher
% and Hayes (2008). 
% 
% INPUTS:
%       x,y = Nx1 vectors to be correlated (both are permuted without
%       replacement)
%
%     Optional Input:
%       'niter' = number of permutations to run. Default = 1000;
%       'flag'   = 'spearmans','pearsons' (default)
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
% Revised and released as BRAVO 2.0 by T. Verstynen (2014)
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.

% Reset random number genrator as a precaution (changing to using Jason's approach)
RandStream('mt19937ar','Seed',sum(100*clock));

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

% Get the initial correlation values of the unpermuted series
switch flag  
 case 'spearmans'
  ctype = 0;
  true_cxy = spear(x(:),y(:));
  
 case 'pearsons'
  ctype = 1;
  r = corrcoef(x(:),y(:));
  true_cxy = r(1,2);
  
 otherwise
  error('Unknown Correlation type');
end;

% Run the bootstrap
boot =  [];
for it = 1:niter
  
  indx = floor(rand(1,n)*(n-1))+1;

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

% Calculate the probability of observing a stronger correlation by chance
if mean(boot) > 0;
    p = length(find(boot < 0))/niter;
else
    p = length(find(boot > 0))/niter;
end;
