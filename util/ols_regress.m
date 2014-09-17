function [betas, resid] = ols_regress(y,x);

% Insert intercept
x = x2fx(x);

% Filter for empty data
ind = find(~isnan(sum([y x],2)));

% Calculates the OLS regression with multiple DVs
betas = inv(x(ind,:)'*x(ind,:))*x(ind,:)'*y(ind,:); 
y_hat = x(ind,:)*betas;
resid = NaN(size(y));
resid(ind,:) = y(ind,:)-y_hat;