function [w_new,X_new,M_new,P_new,D_new,welfare] = func_gen_update(w,X,M,P,D,m)
% T_bar & lambda
T_Bar = m.T.^m.theta./(w.*m.tau);
M_mat = sum(M.*T_Bar,1);
lambda_mat = (M.*T_Bar)./M_mat;

% Generate random J from uniform distribution
rng(123); % set seed
u = sort(rand(m.J, 2));

%initial guess should be CES results
sn = w;
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e6,'MaxIter',1000,'TolFun',1e-10,'TolX',1e-10);
[sn] = fsolve(@(sn)func_sn(m,w,D,sn),sn,options);
sn_star = sn;

c_star = (m.sigma-sn_star.^(m.epsilon/m.sigma))./m.sigma*(m.sigma-1)./m.sigma.*exp((1-sn_star.^(m.epsilon/m.sigma))./m.epsilon).*sn_star.*D.*P;
c_nj = u.*c_star'; % c_n^j = u_j * c_star

s_j = u;
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e6,'MaxIter',1000,'TolFun',1e-10,'TolX',1e-10);
[s_j] = fsolve(@(s_j)func_sj(m,P,D,s_j,c_nj),s_j,options);
s_nj = s_j;

% compute yita ? 虚数,可能是由于initial guess引起了sj是负数，导致这里出现虚数
yita_mat = s_nj.^(m.epsilon/m.sigma).*exp((1-s_nj.^(m.epsilon/m.sigma))./m.epsilon).*s_nj.*c_nj.^(m.theta-1);
yita = sum(yita_mat,1);

% compute nu_P 虚数 same reason
nu_mat = m.theta*m.sigma./(m.sigma-s_nj.^(m.epsilon/m.theta)).*s_nj.*c_nj.^(m.theta);
nu_P = c_star./m.J*sum(nu_mat,1);

% update Dn ? H(x) gammainc function
MT_mat = M.*T_Bar.^m.theta;
H_mat = func_hx(m,s_nj).*c_nj^(m.theta-1);
D_new = D.*(MT_mat.*c_star./m.J*m.theta*sum(H_mat,2))^(1/(1+m.theta));

% update w_i
prod = (1-yita).*lambda_mat.*X;
fmkt_cost = w.*m.F.*c_star.^m.theta*sum(MT_mat,1);
net_prof = yita.*lambda_mat.*X-w.*m.F.*c_star.^m.theta*M.*T_Bar^m.theta;
w_new = prod + fmkt_cost + net_prof;

% update M_i, what is f^e? f^e=1
M_new = net_prof./w./f;

% update P_n
P_new = (P^(-m.theta-1)*nu_P.*MT_mat)^(-1/m.theta);

% update X_n, do we need to update X_n?
X_new = w.*m.L;

% welfare
welfare = w./P;

wnum = sum(w_new);
w_new = w_new/wnum;
X_new = X_new/wnum;
P_new = P_new./wnum;

end
