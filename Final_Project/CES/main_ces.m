
clear
close all
clc

%% Parameters: CES
% general level parameters
m.N = 2;
m.sigma = 4; % substitution rate
m.theta = 4.5; % shape parameter of pareto distribution

% country level parameters
m.T = [2;1]; % scale parateter of pareto distribution
m.L = [15;10]; % labor endowment
m.F = [1;1]; % fixed mkt cost
m.f = [1;1]; % fixed entry cost

%% iceberg trade cost
m.tau = 2*(1-eye(m.N))+eye(m.N);

%% Equiliburium results: CES case
[w,X,M_new,P_new,D] = func_iter(m);

w_ces = w;
M_ces = M_new;
P_ces = P_new;
X_ces = X;
D_ces = D;

disp([w_ces,M_ces,P_ces,X_ces]);