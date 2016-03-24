function Pnii = F2P(F_file, df1, df2, varargin);

% function Pnii = F2P(F_file, df1, df2, varargin);
%

load_type = 'normal'; % Opts: 'normal','untouch'
thresh_type = 'alpha'; 

% Get the variable input parameters
for v=1:2:length(varargin),
    eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
end

Fnii = niiload(F_file, load_type);

F = Fnii.img;

P = 1-fcdf(F_file, df1, df2);

[fp, fn, fe] = fileparts(F_file);

Pnii = Fnii;
Pnii.img = P;
Pnii.hdr.dime.datatype = 16;

ofname = fullfile(fp, ['p_' fn fe]);
niisave(Pnii, ofname, load_type);


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