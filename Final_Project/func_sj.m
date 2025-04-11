function y = func_sj(m,P,D,s_j,c_nj)
y = m.sigma*c_nj./(m.sigma-s_j.^(m.epsilon/m.sigma))-(m.sigma-1)./m.sigma*exp((1-s_j.^(m.epsilon/m.sigma))./m.epsilon).*(D.*P)';
end