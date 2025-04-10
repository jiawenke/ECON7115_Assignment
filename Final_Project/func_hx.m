function y = func_hx(m,s_nj)
a = m.sigma / m.epsilon;
x1 = 1 / m.epsilon;
gamma_1 = gammainc(a,x1,'upper');

x2 = s_nj.^(m.epsilon/m.sigma)/m.epsilon;
gamma_2 = gammainc(a,x2,'upper');


y = 1 + (m.sigma - 1) * exp(1/m.epsilon) * m.epsilon^((m.sigma/m.epsilon)-1) *(gamma_1 - gamma_2);
end
