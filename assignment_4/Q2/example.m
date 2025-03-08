clear
clc
close all

%% Parameters

m.N = 2;
m.barL = 1;
m.sigma = 4;
m.beta = 0.5;
m.alpha = 0.1; % alpha < beta: unique equilibirum

m.barB = ones(m.N,1);
m.barA = [1.2; 1];
m.tau = 1.3*(1-eye(m.N))+eye(m.N);

%% Test s_1

s = -0.2;
[w,L,P,welfare] = func_eqm_iter(s,m);

%% Optimal s_1

ss_vec = (-0.2:0.002:0.09)';
S = length(ss_vec);

wel_vec = ones(S,1);
L_mat = ones(S,m.N);

for ss = 1:S
    s = ss_vec(ss);
    [w,L,P,welfare] = func_eqm_iter(s,m);

    wel_vec(ss) = welfare;
    L_mat(ss,:) = L';

    disp(s)
end

disp('optimal s_1=')
disp(ss_vec(wel_vec==max(wel_vec)))
%%
figure(1)
qx = ss_vec;
h(1) = line(qx, L_mat(:,1), 'Color', 'r', 'LineStyle', '-.', 'LineWidth', 2);
h(2) = line(qx, L_mat(:,2), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);

legend(h, {'$L_1$','$L_2$'},'Location','northwest','Orientation','vertical', 'FontSize',20,'Interpreter','latex')
xlabel('$s_1$', 'FontSize',22, 'Interpreter','latex')
ylabel('Labor', 'FontSize',22)
set(gca,'fontsize',22)
set(gcf, 'Position', get(0, 'Screensize'))
%set(gca,'Xdir','reverse')
print('result/L', '-dpng')

figure(2)
qx = ss_vec;
h = line(qx, wel_vec, 'Color', 'r', 'LineStyle', '-.', 'LineWidth', 2);

xlabel('$s_1$', 'FontSize',22, 'Interpreter','latex')
ylabel('Welfare', 'FontSize',22)
set(gca,'fontsize',22)
set(gcf, 'Position', get(0, 'Screensize'))
%set(gca,'Xdir','reverse')
print('result/welfare', '-dpng')

disp('optimal s_1=')
disp(ss_vec(wel_vec==max(wel_vec)))

%%
function [w_new,L_new,P_new,welfare] = func_eqm_update(w,L,P,s,m)
% Update the equilibrium outcomes by equilibrium conditions
% s = [s_1,s_2,...s_{N-1}]'; s_N is determined by the government budget balance

s_last = -sum(s.*w(1:end-1).*L(1:end-1))/(w(end)*L(end));
s_vec = [s;s_last];

K_mat = (repmat(w./(m.barA.*L.^(m.alpha)),[1,m.N]).*m.tau).^(1-m.sigma); % dim1:i; dim2:n
lambda_mat = K_mat./repmat(sum(K_mat,1),[m.N,1]);

P_new = sum(K_mat,1)'.^(1/(1-m.sigma));

w_new = sum(lambda_mat.*repmat((1+s_vec)'.*w'.*L',[m.N,1]),2)./L; 

B_vec = (m.barB.*(1+s_vec).*w./P).^(1/m.beta);
L_new = m.barL*B_vec/sum(B_vec);

welfare = sum(B_vec).^(m.beta);

wnum = sum(w_new);
w_new = w_new/wnum;
P_new = P_new/wnum;


end

%%
function [w,L,P,welfare] = func_eqm_iter(s,m)
% Solve the equilibrium by iteration
% s = [s_1,s_2,...s_{N-1}]'; s_N is determined by the government budget balance

w = 1/m.N*ones(m.N,1);
L = m.barL/m.N*ones(m.N,1);
P = ones(m.N,1);

diff = 1;
tol = 1e-6;
cc = 0;
ww = 0.1;

while diff>tol && cc < 10000
    [w_new,L_new,P_new,welfare] = func_eqm_update(w,L,P,s,m);
    r = [w;L;P];
    r_new = [w_new;L_new;P_new];

    diff = max(abs(r-r_new));
    cc = cc+1;
    %disp(diff)
    %disp(cc)

    w = ww*w_new + (1-ww)*w;
    L = ww*L_new + (1-ww)*L;
    P = ww*P_new + (1-ww)*P;
end

end


