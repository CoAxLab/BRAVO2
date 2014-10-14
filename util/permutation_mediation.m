function [coeffs, perms] = permutation_mediation(X,Y,M,W,C,varargin)

% function [coeffs perms] = permutation_mediation(X,Y,M,W,C,opts)
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% PERMUTATION_MEDIATION:
% Performs a permutation mediation between two series of data
% to estimate the signficiance of mediator variables M1-Mi on the relationship
% between X & Y. Follows methods reported in Cerin et al. (2006) & Preacher
% and Hayes (2008). 
% 
% BRAVO can now run several types of mediation models:
%
% 1-step Mediation:
%  Y = C'*X + B*M = C'*X + B*A*X
%
% 1-step Moderated Mediation:
%  Y = C'*X + B*M = C'*X + B*(A*X + E*W + F*W*X)
%
% 2-step Mediation:
%  Y = C'*X + B1*M1 + B2*M2  = C'*X + B1*A1*X + B2*(D*M1+A2*X) 
%    = C'*X + B1*A1*X + B2*X*(D*A1+A2) 
%
% 2-step Moderated Mediation (A1 pathway only):
%  Y = C'*X + B1*M1 + B2*M2  = C'*X + B1*(A1*X + E*W + F*W*X) + B2*(D*M1+A2*X) 
%    = C'*X + B1*(A1*X + E*W + F*W*X)  + B2*X*(D*A1+A2) 
%
% Note: For the moderated mediation models BRAVO does not yet implement probing  
%       as specified by Preacher, Rucker, & Hayes (2013). Checke the GitHub
%       repository for updates in the future.
% 
% INPUTS:
%       X,Y,M,W  = independent, dependent, mediator & moderator variables
%       respectively.  X & Y are Nx1 vectors of continous data. M is an 1xP
%       cell array, where P = 1 or 2 depending on how many mediation steps
%       are specified in the model. Each entry in the M cell array should be 
%       an an NxI vector with I = # of mediator variables.
%
%       (Note: Moderated mediation only works on a1-pathways at the moment)
%
%       C = NxL Matrix of covariates (L = # covariates).  If no covariates desired,
%       then give an empty matrix (i.e., [])
%      
%     Optional Input: 
%           'n_iter'     = number of simulations to perform
%           'reg_type'  = type of regression to use: 'ols_regress' (simple OLS, Default)
%                          or 'qr_regress' (QR decomposition)
%
% OUTPUT:  
%       
%       coeffs = Object of path coefficients a, ab, b, c, cprime, & d (moderator) 
%                This is a structure with S fields for each mediation path simulated.
%
%       perms  = Object of simulation arrays of each pathway coefficient from 
%                the permutation tests.
%                This is a structure with p field for each mediation path simulated.
%
% Written by T. Verstynen & A. Weinstein (2011). Modified by T. Verstynen (2013)
% 
% Revised and released as BRAVO 2.0 by T. Verstynen (2014)
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.

% Reset random number genrator as a precaution (changing to using Jason's approach)
RandStream('mt19937ar','Seed',sum(100*clock));

% Define globals used in the subfunction
global reg_type n_paths n_covs n_mediators n_moderators

n_iter = 1000;
reg_type = 'ols_regress';

% Get variable input parameters
for v=1:2:length(varargin),
    eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
end

if ~sum(strcmpi(reg_type,{'ols_regress','qr_regress'}));
    error('Unknown regression type. Options are ols_regress and qr_regress');
end;

% If M & W were sent as matrices, put them into a cell array
if ~iscell(M); M = {M}; end;
if ~iscell(W); W = {W}; end;

% Number of serial indirect paths being modeled
n_paths      = size(M,2);

% Make column vectors
if size(X,2)>size(X,1); X = X'; end;
if size(Y,2)>size(Y,1); Y = Y'; end;
if size(C,2)>size(C,1); C = C'; end;

for p = 1:n_paths
    if size(M{p},2)>size(M{p},1); M{p} = M{p}'; end;
    if ~isempty(W{p}); 
        if size(W{p},2)>size(W{p},1); W{p} = W{p}'; end;
    end;
end;

% Number of working variables
n_covs       = size(C,2);

n_mediators  = cell2mat(cellfun(@size,M,'UniformOutput',0));
n_mediators  = n_mediators(2:2:end);

n_moderators = cell2mat(cellfun(@size,W,'UniformOutput',0));
n_moderators = n_moderators(2:2:end);

warning off;

% Estimate the original coefficients first
coeffs = estimate_pathways(Y,X,M,W,C);

% For summarizing
coeff_list = {'a','b','c_prime','c','ab','d','adb','e','f'};

% Now run the bootstrap
for it = 1:n_iter
  
    % Permute all but the covariates
    nx = X(randperm(length(X)));
    ny = Y(randperm(length(Y)));

    for p = 1:n_paths     
        nm{p} = M{p}(randperm(size(M{p},1)),:);
        if n_moderators(p);
            nw{p} = W{p}(randperm(size(W{p},1)),:);
        else
            nw{p} = [];
        end;
    end;

    % Test the permuated model
    sim_coeffs = estimate_pathways(ny,nx,nm,nw,C);

    %Store all the simulations
    for p = 1:n_paths;
        for c = 1:length(coeff_list);
            eval(sprintf('perms(p).%s(it,:) = sim_coeffs(p).%s;',coeff_list{c},coeff_list{c}));
        end;
    end

end;
warning on



  
