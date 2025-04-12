% clear
% close all
% clc

%% Parameters: CES
% general level parameters
m.N = 2;
m.sigma = 4; % substitution rate
m.theta = 4.5; % shape parameter of pareto distribution
m.epsilon = 2*m.sigma;
% country level parameters
m.T = [2;1]; % scale parateter of pareto distribution
m.L = [15;10]; % labor endowment
m.F = [1;1]; % fixed mkt cost
m.f = [1;1]; % fixed entry cost
m.J = 10; % # of numbers drawn from U
%% iceberg trade cost
m.tau = 2*(1-eye(m.N))+eye(m.N);

%%
[w,X,M_new,P_new,D_new,c_star,welfare,mu_bar] = func_gen_iter(m);
disp([w,M_new,P_new,D_new,c_star,welfare,mu_bar])

% verify c_star
T_Bar = m.T.^(1/m.theta)./(w.*m.tau);
a = 1./T_Bar;
b = (m.sigma-1)/m.sigma*exp(1/m.epsilon).*D_new.*P_new;
if c_star<min(a,b)
    disp('c_star is verified')
else
    disp('c_star is not verified')
end
