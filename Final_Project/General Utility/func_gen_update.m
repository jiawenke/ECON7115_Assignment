function [w_new,X_new,M_new,P_new,D_new,c_star,welfare,mu_bar] = func_gen_update(w,X,M,P,D,m)

%% Generate random J from uniform distribution
rng(123); % set seed
u = sort(rand(m.J,2));

s_n = w;
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e6,'MaxIter',1000,'TolFun',1e-10,'TolX',1e-10);
[s_n] = fsolve(@(s_n)func_sn(m,w,D,X,s_n),s_n,options);
sn_star = s_n;

c_star = ((m.sigma-sn_star.^(m.epsilon/m.sigma))./m.sigma).*((m.sigma-1)./m.sigma).*exp((1-sn_star.^(m.epsilon/m.sigma))./m.epsilon).*D.*P;
c_nj = u.*c_star'; % c_n^j = u_j * c_star

s_j = ones(m.J,2).*sn_star';
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e6,'MaxIter',1000,'TolFun',1e-10,'TolX',1e-10);
[s_j] = fsolve(@(s_j)func_sj(m,P,D,s_j,c_nj),s_j,options);
s_nj = s_j;

% compute yita 
yita_mat_1 = s_nj.^(m.epsilon/m.sigma)./m.sigma.*exp((1-s_nj.^(m.epsilon/m.sigma))./m.epsilon).*s_nj.*c_nj.^(m.theta-1);
yita_mat_2 = exp((1-s_nj.^(m.epsilon/m.sigma))./m.epsilon).*s_nj.*c_nj.^(m.theta-1);
yita = sum(yita_mat_1,1)./sum(yita_mat_2,1);

% compute nu_P
nu_mat = m.theta*m.sigma./(m.sigma-s_nj.^(m.epsilon/m.sigma)).*s_nj.*c_nj.^(m.theta);
nu_P = (c_star'./m.J).*sum(nu_mat,1);

% T_bar & lambda
T_Bar = m.T.^(1/m.theta)./(w.*m.tau);
M_mat = repmat(M,[1,m.N]).*T_Bar.^(m.theta);
lambda_mat = M_mat./repmat(sum(M_mat,1),[m.N,1]);

% update Dn
MT_mat = sum(M_mat,1);
H_mat = func_hx(m,s_nj).*c_nj.^(m.theta-1);
D_new = D.*((MT_mat'.*c_star./m.J*m.theta).*sum(H_mat,1)').^(1/(1+m.theta));

% update w_i
prod_mat = (1-yita).*lambda_mat.*X;
prod = sum(prod_mat,2);
fmkt_cost = w.*m.F.*c_star.^m.theta.*sum(M_mat,1)';
net_prof_mat = yita.*lambda_mat.*X - w.*m.F.*c_star.^m.theta.*M.*T_Bar.^(m.theta);
net_prof = sum(net_prof_mat,2);
disp(net_prof);
w_new = (prod + fmkt_cost + net_prof)./m.L;

% update M_i
M_new = net_prof./w./m.f;

% update P_n
P_new = (P.^(-m.theta-1).*nu_P'.*MT_mat').^(-1/m.theta);

% update X_n
X_new = w.*m.L;

% welfare
welfare = w./P;

% % aggregate markup
mu_mat = s_nj.*c_nj.^m.theta;
mu_bar = sum(mu_mat,1)';

wnum = sum(w_new);
w_new = w_new/wnum;
X_new = X_new/wnum;


end
