function BRAVO_rmANOVA(Y,G,C,mask_file,varargin);

% function BRAVO_mediation(Y,G,C,mask_file,varargin);
% 
% THIS FUNCITON IS IN BETA TESTING; DO NOT USE UNTIL THIS MESSAGE DISAPPEARS

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


out_file  = 'BRAVO_rmANOVA.nii';
load_type = 'untouch'; % Opts: 'normal','untouch'
norm_type = 'zscore';

% Get variable input parameters
for v=1:2:length(varargin),
    eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
end

% Store the parameters in a log object
parameters = struct('outfile',out_file,'load_type',load_type,...
    'norm_type',norm_type);

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

% Assume that Y is a nifti file
nT = size(Y,2);
for t = 1:nT
    eval(sprintf('Y%d=getmat(Y{t},mask_dim,load_type);',t));
end;

% How many unique groups?
groups = unique(G);
groups = groups(find(~isnan(groups))); % Get rid of NaNs
nG = length(groups);

% output storage arrays
store_array = {'meTime','meGroup','intx'};
param_array = {'F','p'};
param_position = [5 6];

for s = 1:length(store_array);
    for p = 1:length(param_array);
        eval(sprintf('%s_%s = NaN(mask_dim);',store_array{s},param_array{p}));
     end;
end;


% Next run the loops
fprintf('\t Running RM ANOVA \n')
old_vox_perc = 0;  % The counter variable

for i = 1:length(good_vox)

    vox_perc = round(100*(i/length(good_vox)));
    % Update the progress if necessary
    if old_vox_perc < vox_perc;
        fprintf(sprintf('  %d',vox_perc));
        old_vox_perc = vox_perc;        
    end;
    
    for t = 1:nT;
        eval(sprintf('tmp = squeeze(Y%d(vx(i),vy(i),vz(i),:));',t));

        if ~isempty(C);
            [~,~,DV(:,t)] = regress(tmp,C);
        else
            DV(:,t) = tmp;
        end;
    end;
 
    vox_Y  = DV(:);
    vox_S  = repmat([1:size(DV,1)],1,nT);
    vox_S  = vox_S(:);
    vox_F1 = repmat([1:nT],size(DV,1),1);
    vox_F1 = vox_F1(:);

    g_array = zeros(size(DV,1),1);
    for g = 1:nG;
        g_array(find(G==groups(g)))=groups(g);
    end;
    vox_F2 = repmat(g_array,nT,1);

    stats = rm_anova2(vox_Y,vox_S,vox_F1,vox_F2,{'Time','Group'});

    for s = 1:length(store_array);
        for p = 1:length(param_array);
            eval(sprintf('%s_%s(vx(i),vy(i),vz(i)) = stats{s+1,%d};',...
                store_array{s},param_array{p},param_position(p)));
        end;
    end;
end;


% Now save everything
for s = 1:length(store_array);
    for p = 1:length(param_array);
        ofname = fullfile(fp,sprintf('%s_%s_%s%s',store_array{s},param_array{p},fn, fe));
        eval(sprintf('%s_%s_nii = mask;',store_array{s},param_array{p}));
        eval(sprintf('%s_%s_nii.img = %s_%s;',store_array{s},param_array{p},store_array{s},param_array{p}));

        eval(sprintf('o_nii = %s_%s_nii;',store_array{s},param_array{p}));

        switch load_type
        case 'normal'
            %eval(sprintf('save_nii(%s_%s_nii,%s);',store_array{s},param_array{p},ofname));
            save_nii(o_nii,ofname);
        case 'untouch'
            %eval(sprintf('save_untouch_nii(%s_%s_nii,%s);',store_array{s},param_array{p},ofname));
            save_untouch_nii(o_nii, ofname);
        otherwise
            error('Unknown saving option');
        end;
    end;
end;

pID_meTime = FDR(meTime_p(good_vox),0.05);
if isempty(pID_meTime); pID_meTime = 0; end;

pID_meGroup = FDR(meGroup_p(good_vox),0.05);
if isempty(pID_meGroup); pID_meGroup = 0; end;

pID_intx = FDR(intx_p(good_vox),0.05);
if isempty(pID_intx); pID_intx = 0; end;

fprintf('\n \n')
fprintf(sprintf('FDR: Main Effect Time = %2.4f \n',pID_meTime));
fprintf(sprintf('FDR: Main Effect Group = %2.4f \n',pID_meGroup));
fprintf(sprintf('FDR: Group x Time = %2.4f \n',pID_intx));

fprintf('\nDone\n')

