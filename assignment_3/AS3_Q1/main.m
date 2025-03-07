clear
clc
close all

%% parameters
x = [0,0];
x1 = x(1);
x2 = x(2);
h = 10^(-15);


%% at point (2,1)
x = [2,1];
x1=x(1);
x2=x(2);
[fwd_x1,fwd_x2,bwd_x1,bwd_x2,ctd_x1,ctd_x2] = func_diff(x1,x2,h);

disp(['Forward diff at point(2,1):',' ',num2str(fwd_x1),' ',num2str(fwd_x2)])
disp(['Backward diff at point(2,1):',' ',num2str(bwd_x1),' ',num2str(bwd_x2)])
disp(['Centered diff at point(2,1):',' ',num2str(ctd_x1),' ',num2str(ctd_x2)])

% automatic differentiation
x = [2,1];
[auto_yx1,auto_yx2] = func_auto(x);
disp(['Automatic diff at point(2,1):',' ',num2str(auto_yx1),' ',num2str(auto_yx2)])

% analytical solution
syms x_1 x_2
y = (x_1^2+x_1^2/x_2-exp(x_2))*log(x_1+exp(x_2));

df1 = diff(y,x_1);
af1_val = double(subs(df1,[x_1,x_2],[x1,x2]));


df2 = diff(y,x_2);
af2_val = double(subs(df2,[x_1,x_2],[x1,x2]));

disp(['Analytical diff at point(2,1):',' ',num2str(af1_val),' ',num2str(af2_val)])


%% at point (5,3)
x = [5,3];
x1=x(1);
x2=x(2);
[fwd_x1,fwd_x2,bwd_x1,bwd_x2,ctd_x1,ctd_x2] = func_diff(x1,x2,h);

disp(['Forward diff at point(5,3):',' ',num2str(fwd_x1),' ',num2str(fwd_x2)])
disp(['Backward diff at point(5,3):',' ',num2str(bwd_x1),' ',num2str(bwd_x2)])
disp(['Centered diff at point(5,3):',' ',num2str(ctd_x1),' ',num2str(ctd_x2)])


% automatic differentiation
x = [5,3];
[auto_yx1,auto_yx2] = func_auto(x);
disp(['Automatic diff at point(5,3):',' ',num2str(auto_yx1),' ',num2str(auto_yx2)])

% analytical solution
   
syms x_1 x_2
y = (x_1^2+x_1^2/x_2-exp(x_2))*log(x_1+exp(x_2));

df1 = diff(y,x_1);
af1_val = double(subs(df1,[x_1,x_2],[x1,x2]));


df2 = diff(y,x_2);
af2_val = double(subs(df2,[x_1,x_2],[x1,x2]));

disp(['Analytical diff at point(5,3):',' ',num2str(af1_val),' ',num2str(af2_val)])

%%
 % Reference: ChatGPT for the code to calculate the analytical solutions