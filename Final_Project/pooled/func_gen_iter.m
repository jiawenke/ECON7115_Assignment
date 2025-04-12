function [w_new,X_new,M_new,P_new,D_new,c_star,welfare,mu_bar] = func_gen_iter(m)
% Solving the equilibrium outcomes by simple iteration
%% initial guess from CES
w = m.w_ces;
M = m.M_ces;
P = m.P_ces;
X = m.X_ces;
D = m.D_ces;
%% iteration
diff = 1;
tol = 1e-6;
ww = 0.5;
cc = 0;

while diff>tol && cc<10000
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
