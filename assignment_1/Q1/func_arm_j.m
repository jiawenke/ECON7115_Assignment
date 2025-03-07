function J = func_arm_j(w,X,tar,m)


kappa = m.tau.*(1+tar);
K_mat = (repmat(w./m.A,[1,m.N]).*kappa).^(1-m.sigma); % dim1:i; dim2:n
P = sum(K_mat,1)'.^(1/(1-m.sigma));
welfare = X./m.L./P;

lambda_mat = K_mat./repmat(sum(K_mat,1),[m.N,1]); % dim1:i; dim2:n

%% Jacobian
J = ones(length(w),length(w));

lambda_new1 = lambda_mat;
sum_ii = zeros(3,3);

sum_ii(1,1)=w(1)*m.L(1)*lambda_new1(1,1)*(1-m.sigma)/w(1)*(1-lambda_new1(1,1));
sum_ii(1,2)=w(2)*m.L(2)*lambda_new1(1,2)*(1-m.sigma)/w(1)*(1-lambda_new1(1,2));
sum_ii(1,3)=w(3)*m.L(3)*lambda_new1(1,3)*(1-m.sigma)/w(1)*(1-lambda_new1(1,3));
sum_ii(2,1)=w(1)*m.L(1)*lambda_new1(2,1)*(1-m.sigma)/w(2)*(1-lambda_new1(2,1));
sum_ii(2,2)=w(2)*m.L(2)*lambda_new1(2,2)*(1-m.sigma)/w(2)*(1-lambda_new1(2,2));
sum_ii(2,3)=w(3)*m.L(3)*lambda_new1(2,3)*(1-m.sigma)/w(2)*(1-lambda_new1(2,3));
sum_ii(3,1)=w(1)*m.L(1)*lambda_new1(3,1)*(1-m.sigma)/w(3)*(1-lambda_new1(3,1));
sum_ii(3,2)=w(2)*m.L(2)*lambda_new1(3,2)*(1-m.sigma)/w(3)*(1-lambda_new1(3,3));
sum_ii(3,3)=w(3)*m.L(3)*lambda_new1(3,3)*(1-m.sigma)/w(3)*(1-lambda_new1(3,3));

J(1,1)=m.L(1)-(sum_ii(1,1)+sum_ii(1,2)+sum_ii(1,3));
J(1,2)=(1-m.sigma)/w(2)*(lambda_new1(1,1)*w(1)*m.L(1)*lambda_new1(2,1)+lambda_new1(1,2)*w(2)*m.L(2)*lambda_new1(2,2)+lambda_new1(1,3)*w(3)*m.L(3)*lambda_new1(2,3));
J(1,3)=(1-m.sigma)/w(3)*(lambda_new1(1,1)*w(1)*m.L(1)*lambda_new1(3,1)+lambda_new1(1,2)*w(2)*m.L(2)*lambda_new1(3,2)+lambda_new1(1,3)*w(3)*m.L(3)*lambda_new1(3,3));
J(2,1)=(1-m.sigma)/w(2)*(lambda_new1(2,1)*w(1)*m.L(1)*lambda_new1(1,1)+lambda_new1(2,2)*w(2)*m.L(2)*lambda_new1(1,2)+lambda_new1(2,3)*w(3)*m.L(3)*lambda_new1(1,3));
J(2,2)=m.L(2)-(sum_ii(2,1)+sum_ii(2,2)+sum_ii(2,3));
J(2,3)=(1-m.sigma)/w(3)*(lambda_new1(2,1)*w(1)*m.L(1)*lambda_new1(3,1)+lambda_new1(2,2)*w(2)*m.L(2)*lambda_new1(3,2)+lambda_new1(2,3)*w(3)*m.L(3)*lambda_new1(3,3));
J(3,1)=(1-m.sigma)/w(1)*(lambda_new1(3,1)*w(1)*m.L(1)*lambda_new1(1,1)+lambda_new1(3,2)*w(2)*m.L(2)*lambda_new1(1,2)+lambda_new1(3,3)*w(3)*m.L(3)*lambda_new1(1,3));
J(3,2)=(1-m.sigma)/w(2)*(lambda_new1(3,1)*w(1)*m.L(1)*lambda_new1(2,1)+lambda_new1(3,2)*w(2)*m.L(2)*lambda_new1(3,2)+lambda_new1(3,3)*w(3)*m.L(3)*lambda_new1(2,3));
J(3,3)=m.L(3)-(sum_ii(3,1)+sum_ii(3,2)+sum_ii(3,3));


end



