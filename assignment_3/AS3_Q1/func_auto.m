function [auto_yx1,auto_yx2] = func_auto(x)

%syms x_1 x_2 y v0 v1 v2 v3 v4 v5 v6 y x_sym;
%x_sym = [x_1,x_2];

%y = (x_1^2+x_1^2/x_2-exp(x_2))*log(x_1+exp(x_2));

v_1 = x(1);
v0 = x(2);
v1 = v_1^2;
v2 = v1/v0;
v3 = exp(v0);
v4 = log(v_1+v3);
v5 = v2+v1-v3;
v6 = v4*v5;
%v6 =y;
auto_yx1 = 1*(v5*(1/(v_1+v3))+v4*(1/v0+2*v_1))*1;
auto_yx2 = 1*(v5*(1/(v_1+v3))*exp(v0)+v4*(1*(-v1/v0^2)-1*exp(v0)))*1;





end
