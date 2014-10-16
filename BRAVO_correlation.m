function BRAVO_correlation(nii_files,regressor,covariates,mask_file,varargin);

% function BRAVO_correlation(nii_files,IV,covariates,mask_file,optlabel,optval);
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% CORRELATION:
% Performs a voxelwise correlation analysis where covariate effects are
% regressed away first and then a bootstrap correlation performed on the 
% relationship between voxel values from the data files and the regressor vector.
% 
% INPUTS:
%       nii_files = cell or character array of input file locations
%       
%       IV        = Column vector (Nx1) of independnent variable
%       
%       covariates = NxC Matrix of covariates.  Default is a row vector if
%                   left empty.
%
%       mask_file  = pointer to a mask file isolating voxels to be analyzed
%                    (voxels > 0).  Must be in same dimensions as input data.
%       
%       Optional Input:
%           method = 'bootstrap' or 'permutation' (Default)
%
%           out_file = Name string for output files (Default 'BRAVO_correlation.nii')
%
%           load_type = Which loader function ('normal','untouched').  Type
%           'help load_nii' and 'help load_untouch_nii' for more info.
%           Default is 'untouch'. NOTE: Using the 'normal option does not
%           always write correctly. Leaving as 'untouch' is highly
%           recommended unless you really know what you are doing.
%
%           corr_type = Which type of correlation to use.  Parametric
%           ('pearsons' Default) or nonparametric ('spearmans')
%
%           n_iter  = Number of permutations to run in the bootstrap
%           (Default 500 iterations).
%
%           norm_type = How to de-mean the data 'zscore'(Default) or
%           'mean0'. Only works on the Y factor
%           
%
% OUTPUT:  Outputs 3 files with postfixes determined by 'out_file' ID.  The
% output file prefix indicates it's value:
%
%       'corr'     = Correlation value (either Pearson's or Spearman's
%       depending on 'corr_type' value).
%
%       'bca_p', 'bca_inv_p' = P-value of bootstrap using a bias corrected
%       and accelerated adjustment. These show the probability that the
%       null mean is larger than the observed effect size (i.e., strong positive
%       effect sizes will have smaller p-values and strong negative
%       effects will have larger p-values). For visualizing in MRICro or similar
%       positive thresholding viewers, use the inv_p file.
%
%       'perc_p', 'perc_inv_p' = P-value of bootstrap using the standard
%       percentile method.  As with teh BCA, the second file is the 1-p values
%       should be used for visualization when thresholding to a particular p-value
%
%       'std'      = Standard deviation of the bootstrap.
% 
% Written by T. Verstynen (2011). Updated 2013
%
% Revised and released as BRAVO 2.0 by T. Verstynen (2014)
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.

method = 'permutation';
out_file  = 'BRAVO_correlation.nii';
load_type = 'untouch'; % Opts: 'normal','untouch'
corr_type = 'pearsons'; % Opts: 'pearsons','spearmans'
n_iter    = 500;
ratio  = 2/3; % With bootstrap method only
norm_type = 'zscore';

% Get variable input parameters
for v=1:2:length(varargin),
    eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
end

% Set the covariates matrix to a column vector of 1's if left blank
if isempty(covariates)
    covariates = ones(size(regressor));
end;

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
        

rOUT = NaN(mask_dim); ppOUT = NaN(mask_dim); sOUT = NaN(mask_dim);
bcap = NaN(mask_dim);

% Next run the loops
fprintf('\t Correlating voxels \n')
old_vox_perc = 0;  % The counter variable

for i = 1:length(good_vox)

    vox_perc = round(100*(i/length(good_vox)));
    % Update the progress if necessary
    if old_vox_perc < vox_perc;
        fprintf(sprintf('  %d',vox_perc));
        old_vox_perc = vox_perc;
    end;
    
    % Extract the voxels data
    series = squeeze(zY(vx(i),vy(i),vz(i),:));
    
    % Regress out the covariates
    [betas, y_resid] = ols_regress(series,covariates);
    
    switch method;
    case 'bootstrap'
        % Run the bootstrap
        [p,r,boot] = bootstrap_correlations(y_resid,regressor,n_iter,corr_type,ratio);
        
        % Calculate percentile significance
        [c_ci, c_p] = ci_percentile(0,boot);
    
        % Calculate BCA significance
        [c_ci_bca, c_p_bca] = ci_bca(0,boot);

    case 'permutation'
        % Run the permuation
        [p,r,boot] = permutation_correlations(y_resid,regressor,n_iter,corr_type);
        
        % Calculate percentile significance
        [c_ci, c_p] = ci_percentile(r,boot);
	c_p = 1-c_p;

        % Calculate BCA significance
        [c_ci_bca, c_p_bca] = ci_bca(r,boot);
	c_p_bca = 1-c_p_bca;

    otherwise
       error(sprintf('Unknown Analysis Method %s',method));

   end;       
 
    % Store data
    rOUT(vx(i),vy(i),vz(i)) = r;
    ppOUT(vx(i),vy(i),vz(i)) = c_p;
    bcapOUT(vx(i),vy(i),vz(i)) = c_p_bca;
    sOUT(vx(i),vy(i),vz(i)) = std(boot);
end;

% Store the new nifti files
rnii = mask; rnii.img = rOUT;
ppnii = mask; ppnii.img = ppOUT;
inv_ppnii = mask; inv_ppnii.img = 1-ppOUT;
bcapnii = mask; bcapnii.img = ppOUT;
inv_bcapnii = mask; inv_bcapnii.img = 1-bcapOUT;
snii = mask; snii.img = sOUT;

% Assign output names
[fp,fn,fe] = fileparts(out_file);
rfile = fullfile(fp,sprintf('corr_%s%s',fn,fe));
ppfile = fullfile(fp,sprintf('perc_p_%s%s',fn,fe));
inv_ppfile = fullfile(fp,sprintf('perc_inv_p_%s%s',fn,fe));
bcapfile = fullfile(fp,sprintf('bca_p_%s%s',fn,fe));
inv_bcapfile = fullfile(fp,sprintf('bcainv_p_%s%s',fn,fe));
sfile = fullfile(fp,sprintf('std_%s%s',fn,fe));

% Save the output accordingly
switch load_type
    case 'normal'
       save_nii(rnii,rfile);
       save_nii(ppnii,ppfile);
       save_nii(inv_ppnii,inv_ppfile);
       save_nii(bcapnii,bcapfile);
       save_nii(inv_bcapnii,inv_bcapfile);
       save_nii(snii,sfile);
    case 'untouch'
       save_untouch_nii(rnii,rfile);
       save_untouch_nii(ppnii,ppfile);
       save_untouch_nii(inv_ppnii,inv_ppfile);
       save_untouch_nii(bcapnii,bcapfile);
       save_untouch_nii(inv_bcapnii,inv_bcapfile);
       save_untouch_nii(snii,sfile); 
    otherwise
        error('Unknown saving option');
end;
    
fprintf('\nDone\n')
return;


