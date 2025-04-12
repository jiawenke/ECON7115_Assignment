
clear
close all
clc

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

%% Equiliburium results: CES case
[w_new,X_new,M_new,P_new,c_star] = func_iter(m);
D = m.sigma/(m.sigma-1);
m.w_ces = w_new;
m.M_ces = M_new;
m.X_ces = X_new;
m.D_ces = D*ones(m.N,1);
m.P_ces = P_new;

disp([m.w_ces,m.M_ces,m.P_ces,m.D_ces]);

%% verify c_star
% pi = w_new.*m.F; % w_n*F_n=cutoff profit
% pi_mat = ((m.sigma*pi)./(P_new.^(m.sigma-1).*X_new)).^(1/(1-m.sigma));
% c_star = ((m.sigma-1)/m.sigma).*(pi_mat);

T_Bar = m.T.^(1/m.theta)./(m.w_ces.*m.tau);
a = T_Bar^(-1);
b = (m.sigma-1)/m.sigma*exp(1/m.epsilon).*m.D_ces.*m.P_ces;
if c_star < min(a,b)
    disp('c_star is verified')
else
    disp('c_star is not verified')
end
%% Non-CES

[w_new,X_new,M_new,P_new,D_new,c_star,welfare,mu_bar] = func_gen_iter(m);
disp([w_new,M_new,P_new,D_new,c_star,welfare,mu_bar])

%% verify c_star
T_Bar = m.T.^(1/m.theta)./(w_new.*m.tau);
a = 1./T_Bar;
b = (m.sigma-1)/m.sigma*exp(1/m.epsilon).*D_new.*P_new;
if c_star<min(a,b)
    disp('c_star is verified')
else
    disp('c_star is not verified')
end
