%% one step GMM
clear
clc
close all

%% Data generation
m.N = 1000; % sample size
m.x = 5*rand(m.N,6); % 6 explanatory variables

theta_true = [1;4.5;3.3;0;-10;-0.1]; % True coefficients
u = 0.1*randn(m.N,1); % error terms

m.y = 1./(1+exp(-m.x*theta_true))+u; % logit
m.z = [m.x, m.x.^2]; % instruments

csvwrite('data.csv',[m.y,m.x,m.z]) % save the data in to csv file

%data: m.x; m.y; m.z
%parameters to estimate: \theta

%% One-step GMM

[~,K]=size(m.z); % K = number of instruments
[~,J] = size(m.x); % J = number of parameters to estimate, thetas
w_mat = eye(K);

theta0 = zeros(J,1); % initial guess

% w0 should be the inverse of the covariance matrix
m0 = func_mm1(theta0,m);
w_mat = (m0'* m0) / size(m0,1);

options = optimoptions('fminunc','Algorithm','quasi-newton','Display','iter',...
    'MaxIterations',500,'MaxFunctionEvaluations',100000,'OptimalityTolerance',1e-8);

[htheta,fval,exitflag,output] = fminunc(@(theta)func_obj_1sgmm(theta,w_mat,m),theta0,options);

disp([theta_true,htheta])


%%
m_final = func_mm1(htheta,m);
S_final = (m_final' * m_final) / size(m_final,1);
W_final = inv(S_final);

theta = htheta;
% Jacobian
[~,J] = size(m.x);

g_3d = ones(K,J,m.N);
for k = 1:K
    for j = 1:J
        for i = 1:m.N
            g_3d(k,j,i) = -exp(-m.x(i,:)*theta)/(1+exp(-m.x(i,:)*theta))^2*m.z(i,k)*m.x(i,j);
        end
    end
end

gy = sum(g_3d,3)/m.N; %dim1:K=#instruments; dim2:J=#parameters

cov_temp = gy'* W_final * gy;
var_cov = inv(cov_temp);

disp(theta)
disp(var_cov)
