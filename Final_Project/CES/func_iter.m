function [w_new,X_new,M_new,P_new,D] = func_iter(m)
% Solving the equilibrium outcomes by simple iteration
%% initial guess
w = ones(m.N,1)/m.N;
X = w;
M = w;
P = w;
%%
diff = 1;
tol = 1e-6;
ww = 0.2;
cc = 0;

while diff>tol && cc<100
    [w_new,X_new,M_new,P_new,D] = func_update(w,X,M,P,m);
    % r = [w;X;M;P]; % why not iterate M and P?
    % r_new = [w_new;X_new;M_new;P_new];
    r = [w;X];
    r_new = [w_new;X_new];
    diff = max(abs(r-r_new));
    cc = cc+1;
    
    w = ww*w_new+(1-ww)*w;
    X = ww*X_new+(1-ww)*X;
    % M = ww*M_new+(1-ww)*M;
    % P = ww*P_new+(1-ww)*P;
end

end
