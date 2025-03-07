clear
clc
close all

%% parameters

m.N=3;
m.tau = 2*(1-eye(m.N))+eye(m.N);
tar = ones(3,3)*0.00;
% tar(1,1)=0;
% tar(2,2)=0;
% tar(3,3)=0;
m.L = [1;2;5];
m.A = [3;1;1];
m.sigma=4;

%% Iteration
% Newton


[w_new,X_new,welfare,lambda_mat,P] = func_eqm_iter(tar,m);

w_a1_newton = w_new;
disp(w_a1_newton)


%% fsolve
w0 = ones(m.N,1)/m.N;
w = w0;
X = w;

options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e6,'MaxIter',1000,'TolFun',1e-10,'TolX',1e-10);
[w,fval] = fsolve(@(w)func_value(w,X,tar,m),w0,options);

disp(w)
disp(fval)
% Normalization
sum = w(1)+w(2)+w(3);
w1 = w(1)/sum;
w2 = w(2)/sum;
w3 = w(3)/sum;
disp(w1)

w_a1_fsolve = [w1;w2;w3];
w_a1_fsolve;
disp(w_a1_fsolve)
disp([w_a1_fsolve,w_a1_newton])

%% If tau_in = 1.2

m.tau = 1.2*(1-eye(m.N))+eye(m.N);

% Use ``Newton``

[w_new,X_new,welfare,lambda_mat,P] = func_eqm_iter(tar,m);

w_b1_newton = w_new;
welfare_a2 = welfare;
disp([welfare_a2]);
