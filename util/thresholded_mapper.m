function outfile = thresholded_mapper(parameter_file, p_file, alpha)

if nargin < 3 || isempty(alpha)
    alpha = 0.025;
end;

coef = load_untouch_nii(parameter_file);
sig  = load_untouch_nii(p_file);


ind = find(sig.img(:)<alpha);

[fp, fn, fe] = fileparts(parameter_file);
outfile = fullfile(fp,sprintf('thresh_%1.3f_%s%s',alpha,fn,fe));

out = coef;
out.img = zeros(size(coef.img));
out.img(ind) = coef.img(ind);

save_untouch_nii(out,outfile);
