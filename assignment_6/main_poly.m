clear
close all
clc

%% Parameters

m.beta = 0.95; % discouting rate
m.alpha = 0.33;
m.gamma = 0.5;
m.xi = 0.2;
m.delta = 0.07;
m.z = 1;

%%
W0 = 1;
m.J = 5; % # of grids
w_vec = (1:m.J)'/m.J*W0;

a = zeros(1,5); %initial guess
diff = 1;
tol = 1e-6;
cc = 0;
ww = 0.5;

while diff>tol&&cc<10000
    v_vec_func = @(w) a(1) + a(2)*(w) + a(3)*(w).^2 + a(4)*(w).^3 + a(5)*(w).^4;
    v_vec = v_vec_func(w_vec);

    v_mat_out = ones(m.J,1);

    g_vec = ones(m.J,1); % policy function: W_{t+1} = g(W_t)
    for j=1:m.J-1
        rhs_vec = ones(j,1);
        for k=1:j % k: tomorrow j: today
            % m.l=1;

            % Get the optimal l*
            n.A = (w_vec(j)^m.alpha*(1-m.alpha))^(-1/m.gamma);
            n.B = w_vec(j);
            n.C = (m.gamma+m.alpha)*(-1/m.gamma);
            n.D = m.delta*w_vec(j)+w_vec(k);

            equation = @(l) n.A * n.B * l^(1 - m.alpha) - l^n.C - n.D * n.A;

            % initial guess
            l_0 = 1;

            % use fsolve
            options = optimoptions('fsolve', 'Display', 'iter');
            m.l = fsolve(equation, l_0, options);
            c_k = (1-m.alpha)*w_vec(j)^m.l-m.delta*w_vec(j)-w_vec(k);
           
            rhs_vec(k) = func_utility_poly(w_vec(j)-w_vec(k),m)+m.beta*v_vec(k);
        end
        v_mat_out(j) = max(rhs_vec);
        g_vec(j) = w_vec(rhs_vec==max(rhs_vec)); % for graph use, here only iterate Value func
    end
    
    
    diff = max(abs(v_vec-v_mat_out));
    cc = cc+1;

    disp(diff)
    disp(cc)
    
    v_vec = ww*v_mat_out + (1-ww)*v_vec;
end

% Fit `a' using OLS
    X = [ones(5, 1), v_vec(1:5), (v_vec(1:5)).^2, (v_vec(1:5)).^3, (v_vec(1:5)).^4]; 
    a = (X' * X) \ (X' * v_mat_out(1:5));
    

disp(a)

%% Stochastic
a = zeros(1,5); %initial guess

z_vec = [0.8;1.2];
pi_mat = [
    0.8, 0.2
    0.4, 0.6]; % 0.8=\pi_LL 0.2=\pi_LH

S = length(z_vec);

diff = 1;
tol = 1e-6;
cc = 0;
ww = 0.5;
%%
while diff>tol&&cc<10000
    v_vec_func = @(w) a(1) + a(2)*(w) + a(3)*(w).^2 + a(4)*(w).^3 + a(5)*(w).^4; % 当前价值函数
    v_mat = repmat(v_vec_func(w_vec(1:5)),[1,S]).*repmat(z_vec',[5,1]);

    v_mat_out = ones(5,S);
    g_mat = ones(5,S); % policy function: W_{t+1} = g(W_t,z_t)
    for j=1:5
        rhs_1 = ones(j,1);
        rhs_2 = ones(j,1);
        for k=1:j

            % Get the optimal l*
            n.A_1 = (z_vec(1) *w_vec(j)^m.alpha*(1-m.alpha))^(-1/m.gamma);
            n.A_2 = (z_vec(2) *w_vec(j)^m.alpha*(1-m.alpha))^(-1/m.gamma);
            n.B_1 = z_vec(1) * w_vec(j);
            n.B_2 = z_vec(2) * w_vec(j);
            n.C = (m.gamma+m.alpha)*(-1/m.gamma);
            n.D = m.delta*w_vec(j)+w_vec(k);
            
            % solve l_* when Z_t=Z_L
            equation = @(l) n.A_1 * n.B_1 * l^(1 - m.alpha) - l^n.C - n.D * n.A;
            l1_0 = 1;
            options = optimoptions('fsolve', 'Display', 'iter');
            m.l_1 = fsolve(equation, l1_0, options);

            % solve l_* when Z_t=Z_H
            equation = @(l) n.A_2 * n.B_2 * l^(1 - m.alpha) - l^n.C - n.D * n.A;
            l2_0 = 1;
            options = optimoptions('fsolve', 'Display', 'iter');
            m.l_2 = fsolve(equation, l2_0, options);
            
            % Bellman equation
            c_k_1 = (1-m.alpha)*w_vec(j)^m.l_1-m.delta*w_vec(j)-w_vec(k);
            c_k_2 = (1-m.alpha)*w_vec(j)^m.l_2-m.delta*w_vec(j)-w_vec(k);
           
            rhs_1(k) = z_vec(1)*func_utility_poly_1(c_k_1,m)+m.beta*(pi_mat(1,1)*v_mat(k,1)+pi_mat(1,2)*v_mat(k,2));
            rhs_2(k) = z_vec(2)*func_utility_poly_2(c_k_2,m)+m.beta*(pi_mat(2,1)*v_mat(k,1)+pi_mat(2,2)*v_mat(k,2));
        end

        v_mat_out(j,1) = max(rhs_1);
        v_mat_out(j,2) = max(rhs_2);
        g_mat(j,1) = w_vec(rhs_1==max(rhs_1));
        g_mat(j,2) = w_vec(rhs_2==max(rhs_2));
    end
    
    diff = max(abs(v_mat(:)-v_mat_out(:)));
    cc = cc+1;
    disp(diff)
    disp(cc)
    
    v_mat = ww*v_mat_out + (1-ww)*v_mat;
end

% Fit `a' using OLS
  X = [ones(5, 1), v_vec(1:5), (v_vec(1:5)).^2, (v_vec(1:5)).^3, (v_vec(1:5)).^4]; 
  a = (X' * X) \ (X' * v_mat_out(1:5));
    
disp(a)

%% Utility functions by z_t

% z_t=1
function y = func_utility_poly(c,m)

y = c.^(1-m.gamma)/(1-m.gamma)-m.l^(1+m.xi)/(1+m.xi);

end

% z_t=z_L
function y = func_utility_poly_1(c,m)

y = c.^(1-m.gamma)/(1-m.gamma)-m.l_1^(1+m.xi)/(1+m.xi);

end

% z_t=z_H
function y = func_utility_poly_2(c,m)

y = c.^(1-m.gamma)/(1-m.gamma)-m.l_2^(1+m.xi)/(1+m.xi);

end
