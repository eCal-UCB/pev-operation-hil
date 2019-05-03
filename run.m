% This script is to run a simulation of an EV charging tariff controller
% that controls charging tariff of two charging options:
% (i) Charging with flexibitliy
% (ii) Charging as soon as possible
%
% EE227C project, May 2019.

%% Initialization
p = init_params();
global par
par = p;

%% Run algorithm -- block coordinate descent
itermax = 1e4;
count = 0;
xk = ones(2*par.N_flex+1,1);            % [soc0, ..., socN, u0, ..., uNm1];
zk = ones(4,1);                         % [z_c, z_uc, y, 1];
vk = 0.1*ones(3,1);                     % [sm_c, sm_uc, sm_y];
par.x0 = xk; par.z0 = zk; par.v0 = vk;
while count < itermax
    % Update x -- fmincon with 
    xk = argmin_x(zk,[],vk);
    
    % Update z -- fmincon
    zk = argmin_z([],xk,vk);
    
    % Update v -- fmincon
    vk = argmin_v(zk,xk,[]);
end

%% Visualization
% TODO: what figures do we want to show?
% TODO: what are the expected results? 