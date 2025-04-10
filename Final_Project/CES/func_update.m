function [w_new,X_new,M_new,P_new,D] = func_update(w,X,M,P,m)

% M_mat = sum(M.*T_Bar^(m.theta),1);
% lambda_mat = (M.*T_Bar)./M_mat;
% X1_mat = lambda_mat.*repmat(X',[m.N,1]); % X_in

T_Bar = m.T.^(1/m.theta)./(w.*m.tau);
M_mat = repmat(M,[1,m.N]).*T_Bar.^(m.theta);
lambda_mat = M_mat./repmat(sum(M_mat,1),[m.N,1]); % dim1:i; dim2:n
X1_mat = lambda_mat.*repmat(X',[m.N,1]); % dim1:i; dim2:n X_in

X2_mat = sum(X1_mat,2); % sum the second dim: n
M_new = (m.sigma-1)/(m.theta*m.sigma) * X2_mat./(w.*m.f); % firm mass

prod = (1-1/m.sigma)*X2_mat; % production wage income
mkt = (m.theta-m.sigma+1)/(m.theta*m.sigma)*X; % mkt fixed cost
net_prof = (m.sigma-1)/(m.theta*m.sigma)*X2_mat; % net profit
w_new = (prod+mkt+net_prof)./m.L;

X_new = w.*m.L; % good market & total expenditure

% try to derive P again
pi = w.*m.F; % w_n*F_n?
pi_mat = (m.sigma*w.*m.F)./(P.^(m.sigma-1).*X);
c_star = (pi_mat).^(1/(1-m.sigma)).*(m.sigma-1)/m.sigma;
P_new = pi.*c_star.^m.theta.*sum(M_mat,1)';

% D is constant, does it need to update?
D = m.sigma/(m.sigma-1);

wnum = sum(w_new);
w_new = w_new/wnum;
X_new = X_new/wnum;
P_new = P_new/wnum;
end
