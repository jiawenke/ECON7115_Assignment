function y = func_hx(m,sn)
a = m.sigma / m.epsilon;
x1 = 1 / m.epsilon;
gamma_1 = gamma(a) * (1 - gammainc(x1, a, 'upper'));

% 计算 Gamma(sigma/epsilon, x^(theta)/epsilon)
x2 = sn.^(m.epsilon) / m.sigma; 
gamma_2 = gamma(a) * (1 - gammainc(x2, a, 'upper'));

% 计算 H(x)
y = 1 + (sigma - 1) * exp(1/epsilon) * epsilon^((sigma-1)/epsilon) * ...
         (gamma_1 - gamma_2);
end