function BRAVO_regression(nii_files,regressor,covariates,contrast,mask_file,varargin);

% function BRAVO_regression(nii_files,IV,covariates,contrast,mask_file,optlabel,optval);
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% REGRESSION:
% Performs a voxelwise regression analysis using bootstrap to estimate the 
% the signficiance of a particular contrast vector (C) of interest.
% 
% INPUTS:
%       nii_files = cell or character array of input file locations
%       
%       IV = Column vector (NxR) of independent variables to regress voxel data
%                   against.
%       
%       covariates = NxC Matrix of covariates.  Default is a row vector if
%                   left empty.
%
%       contrast  = Rx1 contrast vector.  Terms for the covariates and
%                   constant will be added by BRAVO
%
%       mask_file  = pointer to a mask file isolating voxels to be analyzed
%                    (voxels > 0).  Must be in same dimensions as input data.
%       
%       Optional Input: 
%           method = 'bootstrap' or 'permutation' (default)
%
%           out_file = Name string for output files (Default 'BRAVO_regression.nii')
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
%           'mean0'. Only works on the Y factor
%
%           con_type  = What type of contrast to run.  Options are 't' for
%           T-test contrast (i.e., B1-B2 ./ SE(B1-B2)) or 'simple' (Default)
%           contrast that reflects simple combination of regressors (i.e., B1-B2).
%           Future work will allow voxel-wise F-tests.
%           
%           reg_type  = type of regression to use: 'ols_regress' (simple OLS, Default)
%           or 'qr_regress' (QR decomposition)
%
% OUTPUT:  Outputs 3 files with postfixes determined by 'out_file' ID.  The
% output file prefix indicates it's value:
%
%       'con'     = Contrast test value.  If 'con_type' is 't', then this
%       is an image of t-values.  If 'con_type' is 'simple', then this is a
%       difference/sum of regression coefficients (depending on how
%       contrast matrix is setup).
%
%       'bca_p', 'bca_inv_p' = P-value of bootstrap using a bias corrected
%       and accelerated adjustment. These show the probability that the
%       null mean is larger than the observed effect size (i.e., strong positive
%       effect sizes will have smaller p-values and strong negative
%       effects will have larger p-values).  For visualizing in MRICro or similar
%       positive thresholding viewers, use the inv_p file.
%
%       'perc_p', 'perc_inv_p' = P-value of bootstrap using the standard
%       percentile method.  As with teh BCA, the second file is the 1-p values
%       should be used for visualization when thresholding to a particular p-value
%
%       'std'     = Standard deviation of the bootstrap.
%
%       'simulation_permutation_parameters.mat' = A vector of simulation
%       distribution parameters from the permutation or bootstrap
%       (depending on method used).  The order of the vector is [mean_a
%       std_a mean_b std_b mean_ab std_ab mean_c std_c]
%
%      'parameters_log.mat' = An object describing all the parameters
%      used in the analysis run.
% 
% Written by T. Verstynen (2011). Updated in 2013
%
% Revised and released as BRAVO 2.0 by T. Verstynen (2014)
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.

method = 'permutation';
out_file  = 'BRAVO_regression.nii';
load_type = 'untouch'; % Opts: 'normal','untouch'
n_iter    = 500;
norm_type = 'zscore';
con_type  = 't'; % Opts: 't', 'simple' -> just the betas
reg_type  = 'ols_regress'; 

% Get variable input parameters
for v=1:2:length(varargin),
    eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
end

% Store the parameters in a log object
parameters = struct('method',method,'outfile',out_file,'load_type',load_type,...
    'n_iter',n_iter,'norm_type',norm_type,'con_type',con_type,'reg_type',reg_type);

% Get the output path to store the parameters
[fp,fn,fe] = fileparts(out_file);
save(fullfile(fp,'parameters_log.mat'),'parameters');

% Load the mask file to find the voxels of interest
if isempty(mask_file) | ~exist(mask_file,'file')
    error('Please designate a mask file to isolate voxels of interest');
else
    mask = niiload(mask_file,load_type);
    
    % Get the index of all positive voxels
    good_vox = find(mask.img(:)>0);
    
    % Convert to x, y, z
    mask_dim = size(mask.img);
    [vx,vy,vz] = ind2sub(mask_dim,good_vox);
end;

% If no regressor is given
if isempty(contrast)
    contrast = ones(size(regressor,2),1);
    con_type = 'simple';
    fprintf('\t No contrast matrix given.  Defaulting to a vector of ones and averaging betas for analysis.\n')
end;
if size(contrast,2) > size(contrast,1); contrast = contrast'; end;

% Load all the data files
if ischar(nii_files);
    n_files = size(nii_files,1);
elseif iscell(nii_files);
    n_files = length(nii_files);
else
    error('Unknown data file list type.  Has to be an NxP character array or cell array.');
end;

fprintf('\t Loading data files\n')
Y = getmat(nii_files,mask_dim,load_type);

% Normalize the data in the 4th dimension
switch norm_type
    case 'zscore'
        zY = zscore(Y,0,4);
    case 'mean0'
        mY = squeeze(mean(Y,4));
        zY = Y - repmat(mY,[1 1 1 size(Y,4)]);
        
    otherwise
        error('Unknown normalization method');
end;

% Setup outputs
tOUT = NaN(mask_dim); ppOUT = NaN(mask_dim); sOUT = NaN(mask_dim);
bcapOUT = NaN(mask_dim); conOUT = NaN(mask_dim);

% Modify contrast to include the covariates and constant terms
contrast = [contrast; zeros(size(covariates,2),1)];

% Next run the loops
fprintf('\t Regressing voxels \n')
old_vox_perc = 0;  % The counter variable

% Since the design matrix never changes (i.e., brain is always an outcome here)
if isempty(covariates);
    dsX = [regressor];
else
    dsX = [regressor covariates];
end;

for i = 1:length(good_vox)

    vox_perc = round(100*(i/length(good_vox)));
    % Update the progress if necessary
    if old_vox_perc < vox_perc;
        fprintf(sprintf('  %d',vox_perc));
        old_vox_perc = vox_perc;
    end;
    
    % Extract the voxels data
    series = squeeze(zY(vx(i),vy(i),vz(i),:));
    
    % Evaluate the model
    [con, c_p, c_p_bca, m_sim, std_sim] = run_model(dsX, series, contrast, parameters);

    % Store data in matrix form
    conOUT(vx(i),vy(i),vz(i)) = con;
    ppOUT(vx(i),vy(i),vz(i)) = c_p;
    bcapOUT(vx(i),vy(i),vz(i)) = single(c_p_bca);
    sOUT(vx(i),vy(i),vz(i)) = std_sim;
    
    % Store the output
    boot_par(:,i) = [m_sim std_sim];

end;

% Store the new nifti files
connii = mask;      connii.img = conOUT;
ppnii = mask;       ppnii.img = ppOUT;
bcapnii = mask;     bcapnii.img = bcapOUT;
inv_ppnii = mask;   inv_ppnii.img = 1-ppOUT;
inv_bcapnii = mask; inv_bcapnii.img = 1-bcapOUT;
snii = mask;        snii.img = sOUT;

% Assign output names
confile = fullfile(fp,sprintf('con_%s%s',fn,fe));
ppfile = fullfile(fp,sprintf('perc_p_%s%s',fn,fe));
bcapfile = fullfile(fp,sprintf('bca_p_%s%s',fn,fe));
inv_ppfile = fullfile(fp,sprintf('perc_inv_p_%s%s',fn,fe));
inv_bcapfile = fullfile(fp,sprintf('bca_inv_p_%s%s',fn,fe));
sfile = fullfile(fp,sprintf('std_%s%s',fn,fe));

% Save the output accordingly
switch load_type
    case 'normal'
       save_nii(connii,confile);
       save_nii(ppnii,ppfile);
       save_nii(bcapnii,bcapfile);
       save_nii(inv_ppnii,inv_ppfile);
       save_nii(inv_bcapnii,inv_bcapfile);
       save_nii(snii,sfile);
    case 'untouch'
       save_untouch_nii(connii,confile);
       save_untouch_nii(ppnii,ppfile);
       save_untouch_nii(bcapnii,bcapfile);
       save_untouch_nii(inv_ppnii,inv_ppfile);
       save_untouch_nii(inv_bcapnii,inv_bcapfile);
       save_untouch_nii(snii,sfile); 
    otherwise
        error('Unknown saving option');
end;

% Finally save out the bootstrap parameters in the first place
boot_parameters = squeeze(mean(boot_par,2));
save(fullfile(fp,'simulation_distribution_parameters.mat'),'boot_parameters');

fprintf('\nDone\n')
return;


% -----------------------------------------
function [t, c_p, c_p_bca, m_sim, std_sim] = run_model(dsX, series, contrast, params)

switch params.method;
    case 'bootstrap';   
        % Run the bootstrap
        [p, t, ci, sim] = bootstrap_regression(dsX,series,contrast',...
            'stat_type',params.con_type,'n_iter',params.n_iter,...
            'reg_type',params,reg_type);

        % Calculate percentile significance
        [c_ci, c_p] = ci_percentile(0,sim');
    
        % Calculate BCA significance
        [c_ci_bca, c_p_bca] = ci_bca(0,sim');

    case 'permutation'
        % Run the permutation
        [p, t, ci, sim] = permutation_regression(dsX,series,contrast',...
            'stat_type',params.con_type,'n_iter',params.n_iter,...
            'reg_type',params.reg_type);

        % Calculate percentile significance
        [c_ci, c_p] = ci_percentile(t,sim');
        c_p = 1-c_p;

        % Calculate BCA significance
        [c_ci_bca, c_p_bca] = ci_bca(t,sim');
        c_p_bca = 1-c_p_bca;
    
    otherwise
       error(sprintf('Unknown Analysis Method %s',method));
   end;

m_sim = nanmean(sim); std_sim = nanstd(sim);

return;

