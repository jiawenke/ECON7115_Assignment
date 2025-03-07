%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for question a)
function value = func_value(w,X,tar,m)

value = ones(length(w),1);

kappa = m.tau.*(1+tar);
K_mat = (repmat(w./m.A,[1,m.N]).*kappa).^(1-m.sigma); % dim1:i; dim2:n
P = sum(K_mat,1)'.^(1/(1-m.sigma));
welfare = X./m.L./P;

lambda_mat = K_mat./repmat(sum(K_mat,1),[m.N,1]); % dim1:i; dim2:n
Xin_mat = lambda_mat.*repmat(X',[m.N,1]); % dim1:i; dim2:n

value = w.*m.L - sum(Xin_mat./(1+tar),2);

end
