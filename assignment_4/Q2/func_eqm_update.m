%% Equilibruim conditions

function [w_new,X_new,P_new,L_new,welfare] = func_eqm_update(w,X,tarr_vec,L,P,m)


kappa = m.tau.*(1+tarr_vec);
% need repmat. Can't do dot multiply here
K_mat = (w./(m.barA.*L.^m.alpha).*kappa).^(1-m.sigma);
lambda_mat = K_mat./repmat(sum(K_mat,1),[m.R,1]); % dim1:i; dim2:n
Xin_mat = lambda_mat.*repmat(X',[m.R,1]);
P_new = sum(K_mat,1)'.^(1/(1-m.sigma));

% should use sum(numerator,1)
L_sum_1 = L(1)+L(2);
L_sum_2 = L(3)+L(4);
L_new_1 = m.barL.*(m.barB.*X./L./P).^(1./m.mu)/L_sum_1;
L_new_2 = m.barL.*(m.barB.*X./L./P).^(1./m.mu)/L_sum_2;
L_new = [L_new_1(1);L_new_1(2);L_new_2(3);L_new_2(4)];

X_new = w.*L+L_new.*sum(tarr_vec./(1+tarr_vec).*lambda_mat.*sum(Xin_mat,2),1)';

w_new = sum(lambda_mat.*Xin_mat./(1+tarr_vec)',2)./L;

welfare = m.barB.*L.^(-m.mu).*X./L./P;

wnum = sum(w_new);
w_new = w_new/wnum;
P_new = P_new/wnum;
X_new = X_new/wnum;
L_new = L_new/wnum;



end

