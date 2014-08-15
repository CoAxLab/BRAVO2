function  mask = BRAVO_meanmaker(nii_files, varargin);

% function  mask = BRAVO_meanmaker(nii_files, varargin);
%
% BRAVO: Bootstrap Regression Analysis of Voxelwise Observations
%
% MEANMAKER:
% Takes a list of target files and makes a mean image out of them.
%
% INPUTS:
%       nii_files = cell or character array of input file locations
%       
%       Optional Input:
%           load_type = which loader function ('normal','untouched').  Type
%           'help load_nii' and 'help load_untouch_nii' for more info.
%           Default is 'normal'.
%
%           mask_name = name of output mask file. Default is
%           'BRAVO_mean.nii'.
%
% Written by T. Verstnen (2011). Updated 2013.
%
% All code is released under BSD 2-clause license (FreeBSD 9.0).  See
% http://opensource.org/licenses/BSD-2-Clause for more information.

load_type = 'normal'; % Opts: 'normal','untouch'
mask_name = 'BRAVO_mean.nii';

% Get the variable input parameters
for v=1:2:length(varargin),
    eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
end

% Setup the mask output
if iscell(nii_files); mask = niiload(nii_files{1},load_type);
elseif ischar(nii_files);
    mask = niiload(deblank(nii_files(1,:)),load_type);
else
    error('Unknown file pointers');
end;

if size(mask.img,4) > 1; error('Can only accept a list of nifti image files for now'); end;

for f = 1:length(nii_files);
    if ischar(nii_files);
        file = deblank(nii_files(f,:));
    else
        file = nii_files{f};
    end;
     
    % Load the data
    nii = niiload(file,load_type);

    % Store in the data matrix
    if f == 1;
        mask.img = nii.img;
    else
        mask.img = [mask.img + nii.img]./2;
    end;
end;

% save the output
niisave(mask,mask_name,load_type);

return


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

