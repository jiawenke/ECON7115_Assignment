
function y = func_obj_1sgmm(theta,w_mat,m)
% The objective function of the one-step gmm estimator

[hm,~] = func_mm1(theta,m);

y = hm'*w_mat*hm;

end
