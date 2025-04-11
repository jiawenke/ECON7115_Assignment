function [w_new,X_new,M_new,P_new,D] = func_update(w,X,M,P,m)

T_Bar = m.T.^(1/m.theta)./(w.*m.tau);
M_mat = repmat(M,[1,m.N]).*T_Bar.^(m.theta);
lambda_mat = M_mat./repmat(sum(M_mat,1),[m.N,1]); % dim1:i; dim2:n
Xin_mat = lambda_mat.*repmat(X',[m.N,1]); % dim1:i; dim2:n

X2_mat = sum(Xin_mat,2); % Xin sum the second dim: n

% labor market clearing
prod = (1-1/m.sigma)*X2_mat; % production wage income
mkt = (m.theta-m.sigma+1)/(m.theta*m.sigma)*X; % mkt fixed cost
net_prof = (m.sigma-1)/(m.theta*m.sigma)*X2_mat; % net profit
w_new = (prod+mkt+net_prof)./m.L;

M_new = net_prof./(w.*m.f); % firm mass

X_new = w.*m.L; % good market & total expenditure

% Price index P_n
pi = w.*m.F; % w_n*F_n=cutoff profit
pi_mat = (m.sigma*pi)./(P.^(m.sigma-1).*X);
c_star = (pi_mat).^(1/(1-m.sigma)).*(m.sigma-1)/m.sigma;
% P_new = pi.*c_star.^m.theta.*sum(M_mat,1)';
% P_new = m.sigma/(m.sigma-m.theta+1)*(m.sigma/(m.sigma-1))^(1-m.sigma)*c_star.^(m.sigma-m.theta+1).*sum(M_mat,1)';

nu_P = m.theta*(m.sigma/(m.sigma-1))^(1-m.sigma)*P.^m.theta.*1/(m.theta-m.sigma+1).*c_star.^(-m.sigma+m.theta+1);
P_mat = nu_P.*M.*T_Bar.^(m.theta);
P_new = sum(P_mat,1)';

% D is constant
D = m.sigma/(m.sigma-1);

wnum = sum(w_new);
w_new = w_new/wnum;
X_new = X_new/wnum;

end
