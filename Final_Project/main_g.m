
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
m.L = [15;1]; % labor endowment
m.F = [1;1]; % fixed mkt cost
m.f = [1;1];
m.J = 5; % # of numbers drawn from U
%% iceberg trade cost
m.tau = 2*(1-eye(m.N))+eye(m.N);

%%
[w,X,M,P,D] = func_gen_iter(m);
