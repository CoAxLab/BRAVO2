function se = se_calc(resids, x);

if size(x,1)<size(x,2); x = x'; end;

n = size(x,1);
sigma = sqrt( nansum(resids.^2) / (n-2) );

for i = 1:size(x,2);
    se(i) = sqrt( sigma.^2 ./ nansum(x(:,i).^2) );
end;
