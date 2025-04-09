%--------------------main.m------------------------------------------
clear
clc
close all

%% Predetermined arameters
m.alpha = 1/3;
m.delta = 0.06;
m.beta = 0.9;

m.k =(0.01:0.01:10)'; % state space for capital k
m.J = length(m.k); % # of grids

%% Simulated shocks
%{
m.S = 10000;
m.T = 300; % just pick ...?

m.k0_vec = ones(m.S,1); % initial capital stock
m.z0_index = ones(m.S,1); % initial productivity
for s=1:m.S
    m.k0_vec(s) = m.k(randperm(m.J,1));
    m.z0_index(s) = randperm(2,1);
end

m.u_mat = rand(m.S,m.T-1);

save('temp/sim_shock.mat', 'm');
%}

load 'temp/sim_shock.mat'


%% Test
c = 0.1;
piHH = 0.2;
piLL = 0.6;
zgap = 0.3;
theta = [c;piHH;piLL;zgap];

% given any theta, we could simulate m?
[k_mat_sim,d_mat_sim] = func_sim(theta,m);

%-------------------------------------------------------------------------
function [k_mat_sim,d_mat_sim] = func_sim(theta,m)
% Simulate capital and cash flows

m.c = theta(1);
m.piHH = theta(2);
m.piLL = theta(3);
m.zgap = theta(4);

m.zL = 1-m.zgap;
m.zH = 1+m.zgap;

%% Value function iteration

% initial guess for v_H(k) and v_L(k)
vH = m.zH*m.k.^m.alpha;
vL = m.zL*m.k.^m.alpha;

diff = 1;
tol = 1e-6;
cc = 0;
ww = 1;

% get the policy function
while diff>tol && cc<10000
    [vH_out,vL_out,jH,jL] = func_value_update(vH,vL,m);
    r = [vH;vL];
    r_out = [vH_out;vL_out];
    
    diff = max(abs(r-r_out));
    cc = cc+1;
    disp(diff)
    disp(cc)
    
    vH = ww*vH_out + (1-ww)*vH;
    vL = ww*vL_out + (1-ww)*vL;
end

%% Simulations

% state space
z_set = [m.zL;m.zH];

z_mat=ones(m.S,m.T); % simulated productivities: dim1:s;dim2:t
for s=1:m.S
    z_mat(s,1)=z_set(m.z0_index(s));
    for t=2:m.T
        u = m.u_mat(s,t-1);
        if z_mat(s,t-1)==m.zL
            if u <=m.piLL
                z_mat(s,t) = m.zL;
            else
                z_mat(s,t) = m.zH;
            end
        else
            if u <=m.piHH
                z_mat(s,t) = m.zH;
            else
                z_mat(s,t) = m.zL;
            end
        end  
    end
end

k_mat_sim = ones(m.S,m.T);
k_mat_sim(:,1) = m.k0_vec;
d_mat_sim = ones(m.S,m.T-1);

for s=1:m.S
    for t=1:m.T-1
        if z_mat(s,t) == m.zL
            k_mat_sim(s,t+1)=m.k(jL(m.k==k_mat_sim(s,t)));
            i = k_mat_sim(s,t+1)-(1-m.delta)*k_mat_sim(s,t);
            d_mat_sim(s,t)= z_mat(s,t)*k_mat_sim(s,t)^m.alpha-i-m.c/2*i^2/k_mat_sim(s,t);
        else
            k_mat_sim(s,t+1)=m.k(jH(m.k==k_mat_sim(s,t)));
            i = k_mat_sim(s,t+1)-(1-m.delta)*k_mat_sim(s,t);
            d_mat_sim(s,t)= z_mat(s,t)*k_mat_sim(s,t)^m.alpha-i-m.c/2*i^2/k_mat_sim(s,t);
        end
    end
end

end


%-------------------------------------------------------------------------

function [vH_out,vL_out,jH,jL] = func_value_update(vH,vL,m)
% Update the value function
% consume all today
vH_out = ones(m.J,1);
vL_out = ones(m.J,1);
jH = ones(m.J,1); % some records 
jL = ones(m.J,1);

for j=1:m.J
    kk = m.k(m.k>=(1-m.delta)*m.k(j)); % assume I should be positive, non-depreciated stock
    H = length(kk);
    vH_vec = vH(end-(H-1):end);
    vL_vec = vL(end-(H-1):end);
  
    i_vec  = kk-(1-m.delta)*m.k(j); %investment
    dH_vec = m.zH*m.k(j)^m.alpha-i_vec-m.c/2*i_vec.^2/m.k(j);
    dL_vec = m.zL*m.k(j)^m.alpha-i_vec-m.c/2*i_vec.^2/m.k(j);

    % the key parts for the value function iteration
    [vH_out(j),temp_jH] = max(dH_vec+m.beta*(m.piHH*vH_vec+(1-m.piHH)*vL_vec));
    [vL_out(j),temp_jL] = max(dL_vec+m.beta*(m.piLL*vL_vec+(1-m.piLL)*vH_vec));
    
    jH(j) = m.J-H+temp_jH;
    jL(j) = m.J-H+temp_jL;
end

end







