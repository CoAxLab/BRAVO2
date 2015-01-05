function BRAVO_mediation(X,Y,M,covariates,mask_file,varargin);

% function BRAVO_mediation(X,Y,M,covariates,mask_file,optlabel,optval);
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% MEDIATION:
% Performs a voxelwise multiple mediator regression model using bootstrap 
% or permutation tests to estimate the signficiance of mediator variable M 
% on the relationship between X & Y. Follows methods reported in 
% Cerin et al. (2006) & Preacher and Hayes (2008)  & Preacher, Rucker and Hayes (2013)
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
%       X,Y,M, = independent, dependent, mediator and moderator, variables
%       respectively.  These can be either Nx1 column vectors or a pointer
%       to a series of N nifti images.  The pointer can either be an
%       NxP character array of pointers (e.g. spm_select output) or an
%       N-dimensional cell array of pointers to file locations. Dimensons of 
%       ALL images must match the mask file dimensions. For a serial model 
%       M should be a 1x2 cell array with each entry being an NxM array.
%
%       covariates = NxC Matrix of covariates.  If no covariates desired,
%       then give an empty matrix (i.e., [])
%
%       mask_file  = pointer to a mask file isolating voxels to be analyzed
%                    (voxels > 0).  Must be in same dimensions as input data.
%      
%       Optional Input:
%           W  = Moderator arrays for the a-pathways (X->M). The W input should
%           have an identical format as M (i.e., if M is a 1x2 cell array, W is 
%           also a 1x2 cell array). Leave an empty input (i.e., []) if you do not
%           wish to estimate moderating effects on one pathway.
% 
%           method = 'bootstrap' or 'permutation' (default)
%
%           out_file = Name string for output files (Default 'BRAVO_mediation.nii')
%
%           load_type = Which loader function ('normal','untouched').  Type
%           'help load_nii' and 'help load_untouch_nii' for more info.
%           Default is 'untouch'. NOTE: Using the 'normal option does not
%           always write correctly. Leaving as 'untouch' is highly
%           recommended unless you really know what you are doing.
%
%           n_iter  = Number of permutations to run in the bootstrap
%           (Default 500 iterations).
% 
%           norm_type = How to de-mean the data 'zscore'(Default) or
%           'mean0'
%
%           reg_type  = type of regression to use: 'ols_regress' (simple OLS, Default)
%           or 'qr_regress' (QR decomposition)
% 
%
% OUTPUT: BRAVO now outputs several files depending on the model being assessed.  
%     If a serial model is specified, then the mediating path coefficients for 
%     each path will be returned with the prefix "Path1" & "Path2". If only one 
%     indirect pathway step is specified (i.e., a standard X->M-Y indirect 
%     pathway model) then you will not  see the "Path2" prefix. 
%
%     Each output file is a 4 dimensional image where the 4th dimension has
%     as many variables as are specified by the model. For example if you have
%     5 mediating factors in your model, your AB output will be a nifti file
%     with 5 separate images.  The ordering of the images corresponds to the
%     order you specified them in your model. 
%
%     Each pathway coefficient has four output image files.
%
%     1) The map of the statistical coefficients from the model specified. See
%        description above for list of possible coefficients
%
%     2) 'bca_p' = P-value of bootstrap using a bias corrected
%        and accelerated adjustment. This shows the probability that the
%        null mean is larger than the observed effect size (i.e., strong positive
%        effect sizes will have smaller p-values and strong negative
%        effects will have larger p-values).  
%
%     3) 'perc_p', = P-value of bootstrap using the standard
%        percentile method.  NOTE: If you need to use a positive thresholding
%        viewer like MRICroN, use the p_inverter.m function to flip the p-values
%        and threshold your visualization at 1-alpha. For example if you want to
%        find all voxels with p < 0.01, you'd threshold your 'inv_' file, that is
%        output from p_inverter, to 0.99.
%
%     4) '_z' = Z-score of observed regression coefficient against the
%        bootstrap distribution (i.e., [coeff - mean(boot)]/std(boot)). For
%        use in conjunction analysis.
%
%     In addition two matlab files are returned. 
%
%   1) 'simulation_permutation_parameters.mat' = A vector of simulation
%       distribution parameters from the permutation or bootstrap
%       (depending on method used).  The order of the vector is [mean_a
%       std_a mean_b std_b mean_ab std_ab mean_c std_c]
%
%   2) 'parameters_log.mat' = An object describing all the parameters
%      used in the analysis run.
% 
% Written by T. Verstynen & A. Weinstein (2011). Updated by T. Verstynen
% 2012, 2013.  BRAVO 2.0 released 2014
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.


method = 'permutation';
out_file  = 'BRAVO_multimediation.nii';
load_type = 'untouch'; % Opts: 'normal','untouch'
n_iter    = 500;
norm_type = 'zscore';
reg_type  = 'ols_regress'; 
W         = [];
has_mod   = 0;

% Get variable input parameters
for v=1:2:length(varargin),
    eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
end

% Store the parameters in a log object
parameters = struct('method',method,'outfile',out_file,'load_type',load_type,...
    'n_iter',n_iter,'norm_type',norm_type,'reg_type',reg_type);

% Get the output path to store the parameters
[fp,fn,fe] = fileparts(out_file);
save(fullfile(fp,'parameters_log.mat'),'parameters');

% Load the mask file to find the voxels of interest
if isempty(mask_file) || ~exist(mask_file,'file')
    error('Please designate a mask file to isolate voxels of interest');
else
    mask = niiload(mask_file,load_type);
    
    % Get the index of all positive voxels
    good_vox = find(mask.img(:)>0);
    
    % Convert to x, y, z
    mask_dim = size(mask.img);
    [vx,vy,vz] = ind2sub(mask_dim,good_vox);
end;

fprintf('\t Loading data files\n')

% Determine which variables are the nifti data
isX = 0; isY = 0; isM = 0;
if isstr(X) | iscell(X);
    if iscell(X); if isstr(X{1}); isX = 1; end; end;
    if isstr(X); isX=1; end;
end;
if isstr(Y) | iscell(Y);
    if iscell(Y); if isstr(Y{1}); isY = 1; end; end;
    if isstr(Y); isY=1; end;
end;
if isstr(M) | iscell(M);
    if iscell(M); if isstr(M{1}); isM = 1; end; end;
    if isstr(M); isM=1; end;
end;

if isX; X = getmat(X,mask_dim,load_type); end;
if isY; Y = getmat(Y,mask_dim,load_type); end;
if isM; M = getmat(M,mask_dim,load_type); end;

if isM;
    parameters.n_med = 1;
    parameters.n_path = 1;
else
    if ~iscell(M);
        parameters.n_med = size(M,2);
        parameters.n_path = 1;
    else
        parameters.n_path = size(M,2);

        for p = 1:parameters.n_path
            parameters.n_med(p) = size(M{p},2);
        end;
    end;
end;

% Right now set the moderator on the a-path to zero if
% no inputs are given
if isempty (W);
    WV = cell(1,parameters.n_path);
else
    WV = W;
    has_mod = 1;
end;

% Store array
store_array = {'ab','a','b'};
secondary_array = {'d','adb'};


% Setup the full output objects first
for p = 1:parameters.n_path
    for s = 1:length(store_array);
        store_str = store_array{s};
        eval(sprintf('Path%d_%s_perc_p = NaN([mask_dim %d]);',p,upper(store_str),parameters.n_med(p)));
        eval(sprintf('Path%d_%s_bca_p = NaN([mask_dim %d]);',p,upper(store_str),parameters.n_med(p)));
        eval(sprintf('Path%d_%s = NaN([mask_dim %d]);',p,upper(store_str),parameters.n_med(p)));
        eval(sprintf('Path%d_%s_z = NaN([mask_dim %d]);',p,upper(store_str), parameters.n_med(p)));
    end

    switch p
    case 1
        ext_array = {'c_prime'};
    case 2
        ext_array = {'d','adb'};
    end
    for e = 1:length(ext_array);
        store_str = ext_array{e};

        eval(sprintf('%s_perc_p = NaN([mask_dim %d]);',upper(store_str),parameters.n_med(p)));
        eval(sprintf('%s_bca_p = NaN([mask_dim %d]);',upper(store_str),parameters.n_med(p)));
        eval(sprintf('%s = NaN([mask_dim %d]);',upper(store_str),parameters.n_med(p)));
        eval(sprintf('%s_z= NaN([mask_dim %d]);',upper(store_str),parameters.n_med(p)));
    end

    if has_mod
        mod_array = {'f','e'};
        for m = 1:length(mod_array);
            store_str = mod_array{m};
            eval(sprintf('Path%d_%s_perc_p = NaN([mask_dim %d]);',p,upper(store_str),parameters.n_med(p)));
            eval(sprintf('Path%d_%s_bca_p = NaN([mask_dim %d]);',p,upper(store_str),parameters.n_med(p)));
            eval(sprintf('Path%d_%s = NaN([mask_dim %d]);',p,upper(store_str),parameters.n_med(p)));
            eval(sprintf('Path%d_%s_z = NaN([mask_dim %d]);',p,upper(store_str),parameters.n_med(p)));
        end;
    end;
end;


% Next run the loops
fprintf('\t Mediating voxels \n')
old_vox_perc = 0;  % The counter variable

for i = 1:length(good_vox)

    vox_perc = round(100*(i/length(good_vox)));
    % Update the progress if necessary
    if old_vox_perc < vox_perc;
        fprintf(sprintf('  %d',vox_perc));
        old_vox_perc = vox_perc;
    end;
    
    % Extract the voxels data as needed
    if isX; IV = squeeze(X(vx(i),vy(i),vz(i),:));
    else IV = X;
    end;
    
    if isY; DV = squeeze(Y(vx(i),vy(i),vz(i),:));
    else DV = Y;
    end;
    
    if isM; MV = squeeze(M(vx(i),vy(i),vz(i),:));
    else MV = M;
    end;
    
    % Normalize as needed
    IV = norm_data(IV, norm_type);
    DV = norm_data(DV, norm_type);
    MV = norm_data(MV, norm_type);

    % Run the mediation
    [coeffs, tests, sims] = run_model(IV,DV,MV,WV,covariates,parameters);

    % Store data
    for p = 1:parameters.n_path
        for s = 1:length(store_array);
            store_str = store_array{s};
            eval(sprintf('Path%d_%s_perc_p(vx(i),vy(i),vz(i),:) = tests(p).%s_p_perc;',p,upper(store_str),store_str));
            eval(sprintf('Path%d_%s_bca_p(vx(i),vy(i),vz(i),:) = tests(p).%s_p_bca;',p,upper(store_str),store_str));
            eval(sprintf('Path%d_%s(vx(i),vy(i),vz(i),:) = coeffs(p).%s;',p,upper(store_str),store_str));
            eval(sprintf('Path%d_%s_z(vx(i),vy(i),vz(i),:) = tests(p).%s_z;',p,upper(store_str),store_str));
        end

        switch p
        case 1
            ext_array = {'c_prime'};
        case 2
            ext_array = {'d','adb'};
        end
        for e = 1:length(ext_array);
            store_str = ext_array{e};

            eval(sprintf('%s_perc_p(vx(i),vy(i),vz(i),:) = tests(p).%s_p_perc;',upper(store_str),store_str));
            eval(sprintf('%s_bca_p(vx(i),vy(i),vz(i),:) = tests(p).%s_p_bca;',upper(store_str),store_str));
            eval(sprintf('%s(vx(i),vy(i),vz(i),:) = coeffs(p).%s;',upper(store_str),store_str));
            eval(sprintf('%s_z(vx(i),vy(i),vz(i),:) = tests(p).%s_z;',upper(store_str),store_str));
        end

        if has_mod
            mod_array = {'f','e'};
            for m = 1:length(mod_array);
                store_str = mod_array{m};
                eval(sprintf('Path%d_%s_perc_p(vx(i),vy(i),vz(i),:) = tests(p).%s_p_perc;',p,upper(store_str),store_str));
                eval(sprintf('Path%d_%s_bca_p(vx(i),vy(i),vz(i),:) = tests(p).%s_p_bca;',p,upper(store_str),store_str));
                eval(sprintf('Path%d_%s(vx(i),vy(i),vz(i),:) = coeffs(p).%s;',p,upper(store_str),store_str));
                eval(sprintf('Path%d_%s_z(vx(i),vy(i),vz(i),:) = tests(p).%s_z;',p,upper(store_str),store_str));
            end;
        end;

    end;
        
    % Store the output
    for f = 1:size(sims,2);
        sim_parameters{f}(i,:) = sims{f};
    end;

end;

% Assign output names: NOTE: ADD MODERATOR COEFF AS WELL
if has_mod
    par_array = {'AB','A','B','E','F','C_PRIME','D','ADB'};
    max_p1 = 6;
else
    par_array = {'AB','A','B','C_PRIME','D','ADB'};
    max_p1 = 4;
end;

val_array = {'_perc_p','_bca_p','_z',''};

for p = 1:parameters.n_path;
    if p == 1;
        n_params = max_p1;
    else
        n_params = length(par_array);
    end;

    for pr = 1:n_params;
    for v = 1:length(val_array);
        if pr < max_p1
            eval(sprintf('trg = Path%d_%s%s;',p,par_array{pr},val_array{v}));
            str = sprintf('Path%d_%s%s_%s',p,par_array{pr},val_array{v},fn);
        else
            eval(sprintf('trg = %s%s;',par_array{pr},val_array{v}));
            str = sprintf('%s%s_%s',par_array{pr},val_array{v},fn);
        end;
        niisave(mask,trg,str,fp,load_type);
    end;
    end;
end;

% Finally save out the bootstrap parameters in the first place
sim_parameters = cellfun(@mean,sim_parameters,'UniformOutput',0);
save(fullfile(fp,'simulation_distribution_parameters.mat'),'sim_parameters');

fprintf('\nDone\n')

return;

% -----------------------------------------
function [coeffs, hyp_tests, sim_par] = run_model(IV,DV,MV,WV,CV,params);

prime_var = {'a','b','ab','c_prime','e','f'};
secondary_var = {'d','adb'};

switch params.method;
    case 'bootstrap';
    % Run the bootstrap
    [coeffs, sim] = bootstrap_mediation(IV,DV,MV,WV,CV,...
        'n_iter', params.n_iter, 'reg_type', params.reg_type);

    for p = 1:params.n_path;
        for s = 1:length(prime_var)
            eval(sprintf('[~,hyp_tests(p).%s_p_perc] = ci_percentile(0,sim(p).%s);',prime_var{s},prime_var{s}));
            eval(sprintf('[~,hyp_tests(p).%s_p_bca] = ci_bca(0,sim(p).%s);',prime_var{s},prime_var{s}));
            eval(sprintf('hyp_tests(p).%s_z = [nanmean(sim(p).%s)] ./ nanstd(sim(p).%s);',prime_var{s},prime_var{s},prime_var{s}));
        end;

        if p == 2;
        for s = 1:length(secondary_var)
            eval(sprintf('[~,hyp_tests(p).%s_p_perc] = ci_percentile(0,sim(p).%s);',secondary_var{s},secondary_var{s}));
            eval(sprintf('[~,hyp_tests(p).%s_p_bca] = ci_bca(0,sim(p).%s);',secondary_var{s},secondary_var{s}));
            eval(sprintf('hyp_tests(p).%s_z = [nanmean(sim(p).%s)] ./ nanstd(sim(p).%s);',secondary_var{s},secondary_var{s},secondary_var{s}));
        end;
        end;
    end;
 
    case 'permutation';

    % Run the bootstrap
    [coeffs, sim] = permutation_mediation(IV,DV,MV,WV,CV,...
		'n_iter', params.n_iter, 'reg_type', params.reg_type);

    for p = 1:params.n_path;
        for s = 1:length(prime_var)
            eval(sprintf('[~,hyp_tests(p).%s_p_perc] = ci_percentile(coeffs.%s,sim(p).%s);',prime_var{s},prime_var{s},prime_var{s}));
            eval(sprintf('[~,hyp_tests(p).%s_p_bca] = ci_bca(coeffs.%s,sim(p).%s);',prime_var{s},prime_var{s},prime_var{s}));
            eval(sprintf('hyp_tests(p).%s_z = [nanmean(sim.%s) - nanstd(sim.%s)] ./ nanstd(sim.%s);',prime_var{s},prime_var{s},prime_var{s},prime_var{s}));
        end;

        if p == 2;
        for s = 1:length(secondary_var)
            eval(sprintf('[~,hyp_tests(p).%s_p_perc] = ci_percentile(coeffs.%s,sim(p).%s);',secondary_var{s},secondary_var{s},secondary_var{s}));
            eval(sprintf('[~,hyp_tests(p).%s_p_bca] = ci_bca(coeffs.%s,sim(p).%s);',secondary_var{s},secondary_var{s},secondary_var{s}));
            eval(sprintf('hyp_tests(p).%s_z = [nanmean(sim.%s) - nanmean(sim.%s)] ./ nanstd(sim.%s);',secondary_var{s},secondary_var{s},secondary_var{s},secondary_var{s}));
        end;
        end;
    end;

    otherwise
        error(sprintf('Unknown Analysis Method %s',method));
end;


for p = 1:params.n_path
    sim_par{p} = [mean(sim(p).ab) std(sim(p).a) mean(sim(p).b) std(sim(p).c_prime) std(sim(p).d) std(sim(p).adb)];
end;

return;



% -----------------------------------------
function out = norm_data(in,type);

cell_data = 0;
if iscell(in); 
    cell_data = 1;
    cell_size = cellfun(@size,in,'UniformOutput',0);

    in = cell2mat(in);
end;

switch type
    case 'zscore'
        out = [in - repmat(nanmean(in),size(in,1),1)] ./ repmat(nanstd(in),size(in,1),1);
    case 'mean0'
        out = [in - repmat(nanmean(in),size(in,1),1)];
    otherwise
        error('Unknown normalization method');
end;

if cell_data;
    tmp = out; clear out;
    for sz = 1:length(cell_size);
        n_col = cell_size{sz}(2);
        out{sz} = tmp(:,1:n_col);
        tmp = tmp(:,n_col+1:end);
    end;
end;


return;


