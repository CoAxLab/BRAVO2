function mat_data = getmat(nii_files,mask_dim,load_type);

% function mat_data = getmat(nii_files,mask_dim,load_type);
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% GET_MAT:
% General transformation function for NIFTI data into a matlab array
%
% This is a core utility function for all main BRAVO routine.
% 
% Released as BRAVO 2.0 by T. Verstynen (2014)
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.

 
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
