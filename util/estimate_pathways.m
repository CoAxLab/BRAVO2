function coeffs = estimate_pathways(Y,X,M,W,C)

% function coeffs = estimate_pathways(Y,X,M,W,C)
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% ESTIMATE_PATHWAYS:
% Estimate the full regression model pathways for all inputs provide
% including controlling for covariates. The moderator input (W) is
% Only estimated on the X->M (i.e., A) pathways (see Hayes Model #7)
%
% This is a core utility function for the mediation functions.
% 
% Released as BRAVO 2.0 by T. Verstynen (2014)
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.


% Grab globals from parent function
global reg_type n_paths n_covs n_mediators n_moderators

% estimate normal X->Y 
% store in the first path list (if doing a serial model)
[ccoef, r] = ols_regress(Y,[X C]);
coeffs.c  = ccoef(2);

% Empty array for adding in previous mediator paths in the
% serial mediator model
tmp_path = [];
for p = 1:n_paths
    
    % If moderators are included do a separate calculation
    % for the a-pathways
    if ~n_moderators(p)
        design = [X tmp_path C];
        eval(sprintf('[acoef, r] = %s(M{p},design);',reg_type));
     else
        design = [X W{p} repmat(X,1,n_moderators(p)).*W{p} tmp_path C];
        eval(sprintf('[acoef, r] = %s(M{p}, design);',reg_type));
    end;

    % Estimate the b-pathways
    design = [X M{p} tmp_path C];
    eval(sprintf('[primes_coef,r] = %s(Y,design);',reg_type));

    % Store the coefficients 
    coeffs(p).a = acoef(2,:);
    coeffs(p).c_prime = primes_coef(2);
    coeffs(p).b = primes_coef(3:n_mediators(p)+2,:)';
    coeffs(p).ab = coeffs(p).a .* coeffs(p).b;

    % If you're running a serial model store the interaction term 
    % between mediators
    if p > 1;
    	% Replace the old b in the earlier path with the new one 
    	back_b = primes_coef(n_mediators(p)+3:end-n_covs,:);
        coeffs(p-1).b = back_b;

        coeffs(p).d = diag(acoef(3:n_mediators(p)+2,:))';
        coeffs(p).adb = coeffs(p-1).a .* (serial_effect * coeffs(p).d) .* coeffs(p).b; 
        serial_effect = coeffs(p).adb;
        coeffs(p).c = NaN;
    else
    	coeffs(p).d = NaN;
    	coeffs(p).adb = NaN;
    	serial_effect = 1;
    end;


    % Temporary moderator array
    mods = acoef(n_moderators(p)+1+n_mediators(p)+p:end-n_covs,:);
    if ~isempty(mods)
        coeffs(p).e = mods(1,:);
        coeffs(p).f = mods(2,:);
    else
        coeffs(p).e = NaN;
        coeffs(p).f = NaN;
    end;

    % Store the current path to include in the next iteration
    tmp_path = [tmp_path M{p}];
    clear mods design;

end;



