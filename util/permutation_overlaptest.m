function [overlap, boot] = bootstrap_overlaptest(set1,set2,fullset,varargin);

niter = 1000;

% Check for appropriate input;
if nargin < 2
    error('Not enough inputs given');
end;

% Get variable input parameters
for v=1:2:length(varargin),
    eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
end

overlap = length(intersect(set1,set2));

n_set1 = length(set1);
n_set2 = length(set2);

for iter = 1:niter
    
    scram = fullset(randperm(length(fullset)));
    scram_set1 = scram(1:n_set1);
    
    scram = fullset(randperm(length(fullset)));
    scram_set2 = scram(1:n_set2);
    
    boot(iter) = length(intersect(scram_set1,scram_set2));
end;
    
    
