function y = func_sn(m,w,D,X,sn)
sn = w;
y = sn.^(m.epsilon/m.sigma)/m.sigma*(m.sigma-1)/m.sigma.*exp((1-sn.^(m.epsilon/m.sigma))./m.epsilon).*sn.*D.*X-w.*m.F;
end