function [w,X,M_new,P_new,D_new,c_star,welfare,mu_bar] = func_gen_iter(m)
% Solving the equilibrium outcomes by simple iteration
%% initial guess from CES
w = [0.5239;0.4761];
M = [2.5001;1.1666];
P = [0.4688;0.3066];
X = [7.8587;4.7609];
D = [1.3333;1.3333];
%% iteration
diff = 1;
tol = 1e-6;
ww = 0.2;
cc = 0;

while diff>tol && cc<100
    [w_new,X_new,M_new,P_new,D_new,c_star,welfare,mu_bar] = func_gen_update(w,X,M,P,D,m);
    r = [w;X];
    r_new = [w_new;X_new];
    diff = max(abs(r-r_new));
    cc = cc+1;
    %disp(diff)
    %disp(cc)
    
    w = ww*w_new+(1-ww)*w;
    X = ww*X_new+(1-ww)*X;
end

end
