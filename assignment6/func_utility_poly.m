function y = func_utility_poly(c,m)
% Utility function: how to deal with labor?l_t

y = c.^(1-m.gamma)/(1-m.gamma)-m.l^(1+m.xi)/(1+m.xi);

end

% function y = func_poly(c,m)
%%% Utility function
% 
% y = c.^(1-m.gamma)/(1-m.gamma)-l^(1+m.xi)/(1+m.xi);
% 
% end