function [w,X,P,L,welfare] = func_eqm_iter(tarr_vec,m)
% Solving the equilibrium outcomes by simple iteration

% initial guess
w = ones(m.R,1)/m.R;
X = w;
L = w;
P = w;
diff = 1;
tol = 1e-6;
ww = 0.2;
cc = 0;
%%
while diff>tol && cc<1000

   [w_new,X_new,P_new,L_new,welfare] = func_eqm_update(w,X,tarr_vec,L,P,m);

    r = [w;X;L];
    r_new = [w_new;X_new;L_new];
    diff = max(abs(r-r_new));
    cc = cc+1;
    
    w = ww*w_new+(1-ww)*w;
    X = ww*X_new+(1-ww)*X;
    L = ww*L_new+(1-ww)*L;
    P = ww*P_new + (1-ww)*P;
end

end
