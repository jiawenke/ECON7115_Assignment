function [w,X,M,P,D,welfare] = func_gen_iter(m)
% Solving the equilibrium outcomes by simple iteration

w = ones(m.N,1)/m.N;
X = w;
M = w;
P = w;
D = w;
% sn = w;
diff = 1;
tol = 1e-6;
ww = 0.2;
cc = 0;

while diff>tol && cc<100
    [w_new,X_new,M_new,P_new,D_new,welfare] = func_gen_update(w,X,M,P,D,m);
    r = [w;X;M,P,D];
    r_new = [w_new;X_new;M_new,P_new,D_new];
    diff = max(abs(r-r_new));
    cc = cc+1;
    %disp(diff)
    %disp(cc)
    
    w = ww*w_new+(1-ww)*w;
    X = ww*X_new+(1-ww)*X;
end

end
