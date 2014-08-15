function  mask = BRAVO_maskmaker(nii_file, threshold, varargin);

% function  mask = BRAVO_maskmaker(nii_file, threshold, varargin);
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% MASKMAKER:
% Takes a target image file and generates a mask output image.  The mask 
% is binary and reflects a thresholded version of the input
%
% INPUTS:
%       nii_file = cell or character array of input file location
%       
%       threshold = Value to threshold the data with (see Opts)
%
%       Optional Input:
%           load_type = which loader function ('normal','untouched').  Type
%           'help load_nii' and 'help load_untouch_nii' for more info.
%           Default is 'normal'.
%
%           mask_name = name of output mask file. Default is
%           'BRAVO_mask.nii'.
%
%           thresh_type = How to do the thresholding. Greater than
%           threshold value ('gt', default) or less than ('lt').
%
% Written by T. Verstnen (2011). Updated 2013.
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.


load_type = 'normal'; % Opts: 'normal','untouch'
mask_name = 'BRAVO_mask.nii';
thresh_type = 'gt'; % Opts: 'gt'(greater than), 'lt' (less than)

% Get the variable input parameters
for v=1:2:length(varargin),
    eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
end

% load the source data
nii = niiload(nii_file,load_type);

% Prep the mask output
mask = nii;
mask.img = zeros(size(nii.img));

% Get the index of relevant voxels
switch thresh_type
    case 'gt'
        ind = find(nii.img(:) > threshold);
    case 'lt'
        ind = find(nii.img(:) < threshold);
    otherwise
        error('Unknown thresholding option')
end;

% Set to 1s
mask.img(ind) = 1;

% Save the output
niisave(mask, mask_name, load_type);


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
function niisave(nii,name,type);

switch type
    case 'normal'
        save_nii(nii, name);
    case 'untouch'
        save_untouch_nii(nii, name);
    otherwise
        error('Unknown data loader type');
end;
return;
    
