function y = func_sj(m,P,D,s_j,c_nj)
y = m.sigma./(m.sigma-s_j.^(m.epsilon/m.sigma)).*c_nj-(m.sigma-1)./m.sigma*exp((1-s_j.^(m.epsilon/m.sigma))./m.epsilon).*s_j.*(D.*P)';
end