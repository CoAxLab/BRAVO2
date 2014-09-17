function D = BRAVO_sim_fullmodel(varargin)

% Define free variables
N = 150;   % Number of samples
a1 = 0.70; % X  -> M1
a2 = 0.01; % X  -> M2;
b1 = 0.01; % M1 -> Y;
b2 = 0.80; % M2 -> Y;
c  = 0.01; % X  -> Y;
d  = 0.5;  % M1 -> M2;
e  = 0.0;  % moderator on M1
f  = 0.1;  % interaction term on X*W
sigma = 1; % Noise (uniform)
n_sim= 50; % How many simulations to run;

% Get variable input parameters
for v=1:2:length(varargin),
    eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
end

% Define the two "upstream" sources of noise as random variables
w = randn(n_sim,N) + randn(n_sim,N).*sigma;   % Moderator on X->M1
x = randn(n_sim,N) + randn(n_sim,N).*sigma;   % Independent var
m1 = a1*x + e*w + f*(x.*w) + randn(n_sim,N).*sigma; % Mediator 1
m2 = a2*x + d*m1 + randn(n_sim,N).*sigma;           % Mediator 2
y  = c*x  + b1*m1 + b2*m2 + randn(n_sim,N).*sigma;  % Dependent var

% Push to output
D = struct('a1',a1,'a2',a2,'b1',b1,'b2',b2,...
	'c',c,'d',d,'e',e,'f',f,'w',w,...
	'x',x,'y',y,'m1',m1,'m2',m2);
