function out_file = p_inverter(p_file);

% function out_file = p_inverter(p_file);
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% P_INVERTER:
% Simple routine for inverting the p-values for use in 
% visualizers that only threshold in a positive direction
% (e.g., MRICroN, MRICroGL)
% 
% Released as BRAVO 2.0 by T. Verstynen (2014)
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.

% Get file name and set output file
[fp, fn, fe] = fileparts(p_file);
out_file = fullfile(fp,['inv_' fn fe])

% Flip the p value
p = load_untouch_nii(p_file);
i_p = p;
i_p.img = 1-i_p.img;

% Save
save_untouch_nii(i_p, out_file)




