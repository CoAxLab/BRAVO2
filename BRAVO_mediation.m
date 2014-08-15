function BRAVO_mediation(X,Y,M,covariates,mask_file,varargin);

% function BRAVO_mediation(X,Y,M,covariates,mask_file,optlabel,optval);
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% MEDIATION:
% Performs a voxelwise multiple mediator regression model using bootstrap 
% or permutation tests to estimate the signficiance of mediator variable M 
% on the relationship between X & Y. Follows methods reported in 
% Cerin et al. (2006) & Preacher and Hayes (2008). 
% 
% INPUTS:
%       X,Y,M  = independent, dependent and mediator variables
%       respectively.  These can be either Nx1 column vectors or a pointer
%       to a series of N nifti images.  The pointer can either be an
%       NxP character array of pointers (e.g. spm_select output) or an
%       N-dimensional cell array of pointers to file locations. Dimensons of 
%       ALL images must match the mask file dimensions.
%
%       covariates = NxC Matrix of covariates.  If no covariates desired,
%       then give an empty matrix (i.e., [])
%
%       mask_file  = pointer to a mask file isolating voxels to be analyzed
%                    (voxels > 0).  Must be in same dimensions as input data.
%      
%       Optional Input:
%           method = 'bootstrap' or 'permutation' (default)
%
%           out_file = Name string for output files (Default 'BRAVO_mediation.nii')
%
%           load_type = Which loader function ('normal','untouched').  Type
%           'help load_nii' and 'help load_untouch_nii' for more info.
%           Default is 'normal'.
%
%           n_iter  = Number of permutations to run in the bootstrap
%           (Default 500 iterations).
%
%           ratio     = subsampling ratio (Default = 2/3rds). Bootstrap only.
% 
%           norm_type = How to de-mean the data 'zscore'(Default) or
%           'mean0'
%
% OUTPUT:  Outputs 5 files with postfixes determined by out_file ID.  The
% output file prefix indicates it's value:
%
%       'a','b','c','ab'     = Contrast test value for each path in the
%       mediation test.  If multiple mediators are used, the files end
%       with 'm1','m2',...'mi', for the 1-i mediators used.  NOTE: The c 
%       pathway is actually c-prime  
%
%       '_z' = Z-score of observed regression coefficient against the
%       bootstrap distribution (i.e., [coeff - mean(boot)]/std(boot)). For
%       use in conjunction analysis.
%
%       'bca_p', 'bca_inv_p' = P-value of bootstrap using a bias corrected
%       and accelerated adjustment. These show the probability that the
%       null mean is larger than the observed effect size (i.e., strong positive
%       effect sizes will have smaller p-values and strong negative
%       effects will have larger p-values).  For negative effects use the
%       1-p file ('inv_p') file to assess significance.
%
%       'perc_p', 'perc_inv_p' = P-value of bootstrap using the standard
%       percentile method.  As with teh BCA, the second file is the 1-p values
%       should be used when estimating significance on negative effects
%   
%       'simulation_permutation_parameters.mat' = A vector of simulation
%       distribution parameters from the permutation or bootstrap
%       (depending on method used).  The order of the vector is [mean_a
%       std_a mean_b std_b mean_ab std_ab mean_c std_c]
%
%      'parameters_log.mat' = An object describing all the parameters
%      used in the analysis run.
% 
% Written by T. Verstynen & A. Weinstein (2011). Updated by T. Verstynen
% 2012 & 2013
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.


method = 'permutation';
out_file  = 'BRAVO_multimediation.nii';
load_type = 'normal'; % Opts: 'normal','untouch'
n_iter    = 500;
ratio  = 2/3; % With bootstrap method only
norm_type = 'zscore';


% Get variable input parameters
for v=1:2:length(varargin),
    eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
end
% Store the parameters in a log object
parameters = struct('method',method,'outfile',out_file,'load_type',load_type,...
    'n_iter',n_iter,'norm_type',norm_type,'ratio',ratio);

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
% see which variable is the nifti data
isX = 0; isY = 0; isM = 0;
if ~isnumeric(X); X = getmat(X,mask_dim,load_type); isX = 1; end;
if ~isnumeric(Y); Y = getmat(Y,mask_dim,load_type); isY = 1; end;
if ~isnumeric(M); M = getmat(M,mask_dim,load_type); isM = 1; end;

% Setup output variables for storage
if size(M,4)>1; n_meds = 1; else; n_meds = size(M,2); end;

ab_pOUT = NaN([mask_dim n_meds 2]); 
a_pOUT = NaN([mask_dim n_meds 2]); 
b_pOUT = NaN([mask_dim n_meds 2]); 
c_pOUT = NaN([mask_dim n_meds 2]); 
abOUT = NaN([mask_dim n_meds]); aOUT = NaN([mask_dim n_meds]); bOUT = NaN([mask_dim n_meds]); 
ab_zOUT = NaN([mask_dim n_meds]); a_zOUT = NaN([mask_dim n_meds]); b_zOUT = NaN([mask_dim n_meds]); 
cOUT = NaN(mask_dim);  
c_zOUT = NaN(mask_dim);  

% Store array
store_array = {'ab','a','b'};

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

    switch method;
    case 'bootstrap';
        % Run the bootstrap
        [coeffs, ab_boot, a_boot, b_boot, c_boot] = bootstrap_mediation(IV,DV,MV,covariates,n_iter,ratio);

        % Calculate percentile significance
        [ab_ci, ab_p] = ci_percentile(0,ab_boot);
        [a_ci, a_p] = ci_percentile(0,a_boot);
        [b_ci, b_p] = ci_percentile(0,b_boot);
        [c_ci, c_p] = ci_percentile(0,c_boot);
    
        % Calculate BCA significance
        [ab_ci_bca, ab_p_bca] = ci_bca(0,ab_boot);
        [a_ci_bca, a_p_bca] = ci_bca(0,a_boot);
        [b_ci_bca, b_p_bca] = ci_bca(0,b_boot);
        [c_ci_bca, c_p_bca] = ci_bca(0,c_boot);
    
        % Estimate the standard deviations
        ab_z = [nanmean(ab_boot')] ./ nanstd(ab_boot');
        a_z  = [nanmean(a_boot')] ./ nanstd(a_boot');
        b_z =  [nanmean(b_boot')] ./ nanstd(b_boot');
        c_z =  [nanmean(c_boot')] ./ nanstd(c_boot');
    
    case 'permutation';
        % Run the bootstrap
        [coeffs, ab_boot, a_boot, b_boot, c_boot] = permutation_mediation(IV,DV,MV,covariates,n_iter);

        % Calculate percentile significance
        [ab_ci, ab_p] = ci_percentile(coeffs.ab,ab_boot);
        [a_ci, a_p] = ci_percentile(coeffs.a,a_boot);
        [b_ci, b_p] = ci_percentile(coeffs.b,b_boot);
        [c_ci, c_p] = ci_percentile(coeffs.cprime,c_boot);
        ab_p = 1-ab_p; a_p = 1-a_p; b_p = 1-b_p; c_p = 1-c_p;

        % Calculate BCA significance
        [ab_ci_bca, ab_p_bca] = ci_bca(coeffs.ab,ab_boot);
        [a_ci_bca, a_p_bca] = ci_bca(coeffs.a,a_boot);
        [b_ci_bca, b_p_bca] = ci_bca(coeffs.b,b_boot);
        [c_ci_bca, c_p_bca] = ci_bca(coeffs.cprime,c_boot);
        ab_p_bca = 1-ab_p_bca; a_p_bca = 1-a_p_bca; b_p_bca = 1-b_p_bca; c_p_bca = 1-c_p_bca;


        % Estimate the standard deviations
        ab_z = [coeffs.ab' - nanmean(ab_boot')] ./ nanstd(ab_boot');
        a_z  = [coeffs.a' - nanmean(a_boot')] ./ nanstd(a_boot');
        b_z =  [coeffs.b' - nanmean(b_boot')] ./ nanstd(b_boot');
        c_z =  [coeffs.c' - nanmean(c_boot')] ./ nanstd(c_boot');
        
    otherwise
        error(sprintf('Unknown Analysis Method %s',method));
        
    end;

    % Store data
    for s = 1:length(store_array);
        store_str = store_array{s};
        eval(sprintf('%s_pOUT(vx(i),vy(i),vz(i),:,1) = %s_p;',store_str,store_str));
        eval(sprintf('%sOUT(vx(i),vy(i),vz(i),:) = coeffs.%s;',store_str,store_str));
        eval(sprintf('%s_pOUT(vx(i),vy(i),vz(i),:,2) = %s_p_bca;',store_str,store_str));
        eval(sprintf('%s_zOUT(vx(i),vy(i),vz(i),:) = %s_z;',store_str,store_str));
    end;
        
    % Finally save the original direct path
	c_pOUT(vx(i),vy(i),vz(i),1,1) = c_p;
    c_pOUT(vx(i),vy(i),vz(i),1,2) = c_p_bca;
    cOUT(vx(i),vy(i),vz(i)) = coeffs.cprime;

    % Store the output
    boot_par(:,i) = [mean(a_boot) std(a_boot) mean(b_boot) std(b_boot) mean(ab_boot) std(ab_boot) mean(c_boot) std(c_boot)];

end;

% Assign output names
par_array = {'ab','a','b','c'};
val_array = {'_p','_inv_p','_z',''};

% First invert the p-value array for ease of display in viewers like MRICRON
for s = 1:length(par_array);
        par = par_array{s};
        str = sprintf('%s_inv_pOUT = 1-%s_pOUT;',par,par);
        eval(str);
end;

for p = 1:length(par_array);
    for v = 1:length(val_array);
        if v < 3;
            isp = 1;
        else
            isp = 0;
        end;

        eval(sprintf('trg = %s%sOUT;',par_array{p},val_array{v}));
        str = sprintf('%s%s_%s',par_array{p},val_array{v},fn);
        niisave(mask,trg,str,fp,load_type,isp);
    end;
end;

% Finally save out the bootstrap parameters in the first place
boot_parameters = squeeze(mean(boot_par,2));
save(fullfile(fp,'simulation_distribution_parameters.mat'),'boot_parameters');

fprintf('\nDone\n')

return;



% -----------------------------------------
function mat_data = getmat(nii_files,mask_dim,load_type);
 
if ischar(nii_files);
    n_files = size(nii_files,1);
elseif iscell(nii_files);
    n_files = length(nii_files);
else
    error('Unknown data file list type.  Has to be an NxP character array or cell array.');
end;

% If it's a 4-D file, then just do that one;
if n_files == 1;
    
    if ischar(nii_files); file = deblank(nii_files(1,:)); 
    else
        file = nii_files{1};
    end;

    nii = niiload(file,load_type);
    mat_data = nii.img;
else
    % Otherwise loop through and load each file
    % Get the data files setup
    mat_data = NaN([mask_dim n_files]);
    for f = 1:n_files

        if ischar(nii_files);
            file = deblank(nii_files(f,:));
        else
            file = nii_files{f};
        end;

        % Load the data
        nii = niiload(file,load_type);

        % Stop if the image has incorrect dimensions
        if sum(size(nii.img)-mask_dim); error(sprintf('Image %d does not match mask image dimensions',f)); end;

        % Store in the data matrix
        mat_data(:,:,:,f) = nii.img;
    end;
end;

return;


% -----------------------------------------
function out_nii = niiload(file,type);

switch type
    case 'normal'
        out_nii = load_nii(file);
    case 'untouch'
        out_nii = load_untouch_nii(file);
    otherwise
        error('Unknown data loader type');
end;
return;


% -----------------------------------------
function niisave(mask,mat,namestr,path,type,isp);

if nargin < 6 || isempty(isp)
    isp = 0;
end;

n=size(mat,4);

if isp
    for f = 1:n
        
        if n > 1;
            name1 = fullfile(path,['perc_' sprintf([namestr '_m%d'],f) '.nii']);
            name2 = fullfile(path,['bca_' sprintf([namestr '_m%d'],f) '.nii']);

        else
            name1 = fullfile(path,['perc_' namestr '.nii']);
            name2 = fullfile(path,['bca_' namestr '.nii']);
        end;

        mat1 = squeeze(mat(:,:,:,f,1));
        mat2 = squeeze(mat(:,:,:,f,2));
        nii1 = mask;
        nii2 = mask;
        
        % Make sure the data format is appropriate
        nii1.img = single(mat1);
        nii1.hdr.dime.datatype = 16;
        nii1.hdr.dime.bitpix = 32;
        
        % Make sure the data format is appropriate
        nii2.img = single(mat2);
        nii2.hdr.dime.datatype = 16;
        nii2.hdr.dime.bitpix = 32;
        
        switch type
            case 'normal'
                save_nii(nii1, name1);
                save_nii(nii2, name2);
            case 'untouch'
                save_untouch_nii(nii1, name1);
                save_untouch_nii(nii2, name2);
            otherwise
                error('Unknown data loader type');
        end;
    end;
else
    for f = 1:n
        
        if n > 1;
            name = fullfile(path,[sprintf([namestr '_m%d'],f) '.nii']);
        else
            name = fullfile(path,[namestr '.nii']);
        end;
        nmat = squeeze(mat(:,:,:,f));
        nii = mask;
        
        % Make sure the data format is appropriate
        nii.img = single(nmat);
        nii.hdr.dime.datatype = 16;
        nii.hdr.dime.bitpix = 32;
        
        switch type
            case 'normal'
                save_nii(nii, name);
            case 'untouch'
                save_untouch_nii(nii, name);
            otherwise
                error('Unknown data loader type');
        end;

    end;
end;    

return

% -----------------------------------------
function out = norm_data(in,type);

switch type
    case 'zscore'
        out = [in - repmat(nanmean(in),size(in,1),1)] ./ repmat(nanstd(in),size(in,1),1);
    case 'mean0'
        out = [in - repmat(nanmean(in),size(in,1),1)];
    otherwise
        error('Unknown normalization method');
end;




