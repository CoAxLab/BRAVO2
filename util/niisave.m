function niisave(mask,mat,namestr,path,type);

% function niisave(mask,mat,namestr,path,type);
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% NIISAVE:
% General saver function for NIFTI data.
%
% This is a core utility function for all main BRAVO routine.
% 
% Released as BRAVO 2.0 by T. Verstynen (2014)
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.


% Get number of dimensions (1 per mediating variable per path)
n=size(mat,4);

% Get the mask nii object
nii = mask;

if n > 1;
    % Reset the nifti header to be 4-d if needed
    nii.hdr.dime.dim(1) = 4;
    nii.hdr.dime.dim(5) = n;
end;

% define the name
name = fullfile(path,[namestr '.nii']);

% Make sure the data format is appropriate
nii.img = single(mat);
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

return
