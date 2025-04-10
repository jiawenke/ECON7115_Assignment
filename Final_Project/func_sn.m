function y = func_sn(m,w,D,X,s_n)
y = (s_n.^(m.epsilon/m.sigma)./m.sigma)*((m.sigma-1)/m.sigma).*exp((1-s_n.^(m.epsilon/m.sigma))./m.epsilon).*s_n.*D.*X-w.*m.F;
end
