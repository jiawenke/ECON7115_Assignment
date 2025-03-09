clear
close all
clc

%% Parameters

m.beta = 0.95; % discouting rate
m.alpha = 0.33;
m.gamma = 0.5;
m.xi = 0.2;
m.delta = 0.07;
m.l = 1;
m.z = 1;

W0 = 1;
m.J = 200; % # of grids
w_vec = (1:m.J)'/m.J*W0;

%%
a = zeros(1,5);
diff = 1;
tol = 1e-6;
cc = 0;
ww = 0.5;

while diff>tol&&cc<10000
    v_vec_func = @(w) a(1) + a(2)*(w) + a(3)*(w).^2 + a(4)*(w).^3 + a(5)*(w).^4; % polynomial of Value function
    v_vec = v_vec_func(w_vec);

    v_vec_out = ones(m.J,1);
    %g_vec = ones(m.J,1); % policy function: W_{t+1} = g(W_t)
    for j=1:m.J
        rhs_vec = ones(j,1);
        for k=1:j % k: tomorrow j: today
            c_k = w_vec(j)^m.alpha*m.z*m.l^(1-m.alpha)-m.delta*w_vec(j)-w_vec(k);
            rhs_vec(k) = func_utility_poly(c_k,m)+m.beta*v_vec(k);
        end
        v_vec_out(j) = max(rhs_vec);
        %g_vec(j) = w_vec(rhs_vec==max(rhs_vec)); % for graph use, here only iterate Value func
    end
    
    % Fit a using OLS
    X = [ones(m.J, 1), w_vec, w_vec.^2, w_vec.^3, w_vec.^4]; % 4th degree polynomial
    a = (X' * X) \ (X' * v_vec_out); % OLS approximation
    
    diff = max(abs(v_vec-v_vec_out));
    cc = cc+1;

    disp(diff)
    disp(cc)
    
    v_vec = ww*v_vec_out + (1-ww)*v_vec;
end

disp(a)
