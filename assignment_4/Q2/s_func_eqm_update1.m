%% Equilibruim conditions

function [w_new,X_new,P_new,L1_new,welfare] = s_func_eqm_update1(w,X,tar_2,L1,P,m)
tar_1 = 3650;
tarr_vec = zeros(m.R,m.R);
tarr_vec(3,1)=tar_1;
tarr_vec(3,2)=tar_1;
tarr_vec(4,1)=tar_1;
tarr_vec(4,2)=tar_1;

tarr_vec(1,3)=tar_2;
tarr_vec(1,4)=tar_2;
tarr_vec(2,3)=tar_2;
tarr_vec(2,4)=tar_2;

kappa = m.tau.*(1+tarr_vec);
K_mat = (repmat(w./m.barA*(m.barL^m.alpha),[1,m.R]).*kappa).^(1-m.sigma);
lambda_mat = K_mat./repmat(sum(K_mat,1),[m.R,1]); % dim1:i; dim2:n
Xin_mat = lambda_mat.*repmat(X',[m.R,1]);
P_new = sum(K_mat,1)'.^(1/(1-m.sigma));


L1_new = (m.barB*X./L1./P).^(1./m.mu);

X_new = w.*L1+L1_new.*sum(tarr_vec./(1+tarr_vec).*lambda_mat.*sum(Xin_mat,2),1)';

w_new = sum(lambda_mat.*Xin_mat./(1+tarr_vec)',2)./L1;

welfare = sum(L1_new).^(1./m.mu);


wnum = sum(w_new);
w_new = w_new/wnum;
P_new = P_new/wnum;

end

%%
% %% Equilibruim conditions
% 
% function [w_new,X_new,P_new,L1_new,welfare] = s_func_eqm_update1(w,X,tar,L1,P,m)
% 
% %tar_2 = 0;
% %tarr_vec = [tar;tar;tar_2;tar_2];
% tarr_vec = zeros(m.R,m.R);
% tarr_vec(3,1)=tar;
% tarr_vec(3,2)=tar;
% tarr_vec(4,1)=tar;
% tarr_vec(4,2)=tar;
% 
% kappa = m.tau.*repmat((1+tarr_vec),[1,m.R]);
% K_mat = (repmat(w./m.barA*(m.barL^m.alpha),[1,m.R]).*kappa).^(1-m.sigma);
% lambda_mat = K_mat./repmat(sum(K_mat,1),[m.R,1]); % dim1:i; dim2:n
% Xin_mat = lambda_mat.*repmat(X',[m.R,1]);
% %Xin_mat = lambda_mat.*X';
% P_new = sum(K_mat,1)'.^(1/(1-m.sigma));
% 
% 
% % L1 L2 matrix
% L1_mat = (m.barB*X./L1./P).^(1/0.5);
% L1_new = m.barL*L1_mat./repmat(sum(L1_mat,1),[m.R,1]);
% %L1_new = m.barL*L1_mat./sum(L1_mat,1);
% 
% %L2_mat = (m.barB*X./L./P).^(1/m.mu2);
% %L2_new = m.barL*L2_mat./repmat(sum(L2_mat,1),[m.R,1]);
% 
% X_new = w.*L1+L1_new.*sum(tarr_vec./(1+tarr_vec).*lambda_mat.*sum(Xin_mat,2),1)';
% 
% w_new = sum(lambda_mat.*repmat(X'/(1+tarr_vec)',[m.R,1]),2)./L1;
% 
% welfare = sum(L1_mat).^(m.mu);
% 
% 
% wnum = sum(w_new);
% w_new = w_new/wnum;
% P_new = P_new/wnum;
% 
% end
