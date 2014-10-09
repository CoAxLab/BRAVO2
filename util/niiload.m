function out_nii = niiload(file,type);

% function out_nii = niiload(file,type);
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% NIILOAD:
% General loader function for NIFTI data.
%
% This is a core utility function for all main BRAVO routine.
% 
% Released as BRAVO 2.0 by T. Verstynen (2014)
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.


switch type
    case 'normal'
        out_nii = load_nii(file);
    case 'untouch'
        out_nii = load_untouch_nii(file);
    otherwise
        error('Unknown data loader type');
end;

return;
