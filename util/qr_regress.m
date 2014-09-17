function [betas, resid] = qr_regress(y,x);

% Insert intercept
x = x2fx(x);

% Filter for empty data
ind = find(~isnan(sum([y x],2)));

% Calculate the betas using QR decomposition
[Q, R] = qr(x(ind,:),0);
betas = R \ (Q' * y(ind,:));
y_hat = x(ind,:)*betas;
resid = NaN(size(y));
resid(ind,:) = y(ind,:)-y_hat;

