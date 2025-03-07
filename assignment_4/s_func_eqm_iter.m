function [w,X,P,L1,welfare] = s_func_eqm_iter(tar,m)
% Solving the equilibrium outcomes by simple iteration

% initial guess
w = ones(m.R,1)/m.R;
X = w;
L = w;
L1 =w;
P = w;
diff = 1;
tol = 1e-6;
ww = 0.2;
cc = 0;

while diff>tol && cc<100
   [w_new,X_new,P_new,L1_new,welfare] = s_func_eqm_update(w,X,tar,L1,P,m);

    r = [w;X;L];
    r_new = [w_new;X_new;L1_new];
    diff = max(abs(r-r_new));
    cc = cc+1;
    
    w = ww*w_new+(1-ww)*w;
    X = ww*X_new+(1-ww)*X;
    L1 = ww*L1_new+(1-ww)*L1;
    %L2 = ww*L2_new+(1-ww)*L2;
    P = ww*P_new + (1-ww)*P;
end

end
