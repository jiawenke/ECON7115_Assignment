clear
clc
close all

%% Parameters
m.N=3;
m.tau = 2*(1-eye(m.N))+eye(m.N);
tar = ones(3,3)*0.00;
% tar(1,1)=0;
% tar(2,2)=0;
% tar(3,3)=0;
m.L = [1;2;5];
m.A = [3;1;1];
m.sigma=4;

%%
% %% Baseline equilibrium, tarrif=0
% simple iteration
% 
 [w,X,P,welfare,lambda_mat] = func_eqm_iter(tar,m);
% 





%% tarrif=0.05

tar = ones(m.N,m.N)*0.05;
tar(1,1) = 0;
tar(2,2) = 0;
tar(3,3) = 0;

[w,X,P,welfare,lambda_mat] = func_eqm_iter(tar,m);

GE_005 = 1-m.sigma-lambda_mat;
disp(GE_005)


%% tarrif=0.25

tar = ones(m.N,m.N)*0.25;
tar(1,1) = 0;
tar(2,2) = 0;
tar(3,3) = 0;

[w,X,P,welfare,lambda_mat] = func_eqm_iter(tar,m);

GE_025 = 1-m.sigma-lambda_mat;
disp(GE_025)



%}

