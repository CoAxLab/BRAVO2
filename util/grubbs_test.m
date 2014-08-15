function [g, z, g_thresh] = grubbs_test(x,alpha);

if nargin < 2 | isempty(alpha);
    alpha = 0.05;
end;

is_outlier = 1;

g = zeros(size(x));

while is_outlier
    g_dist = abs(x-nanmean(x)) ./ nanstd(x);
    z = max(g_dist);
    n = sum(~isnan(x));
    tsq = tinv(1-(alpha/(2*n)),n-2).^2;
    g_thresh = [(n-1)./sqrt(n)] * sqrt(tsq./(n-2+tsq));

    if z > g_thresh;
        is_outlier = 1;
        g(find(abs(g_dist) == z))=1;
        x(find(abs(g_dist) == z))=NaN;
    else
        is_outlier = 0;
    end;
end;

