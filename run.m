% This script is to run a simulation of an EV charging tariff controller
% that controls charging tariff of two charging options:
% (i) Charging with flexibitliy
% (ii) Charging as soon as possible
%
% EE227C project, May 2019.

%% Initialization
par = init_params();


%% Run algorithm -- block coordinate descent
% TODO: create matlab function for each choice -- to run fmincon
itermax = 1e4;
count = 0;
xk = []; % TODO: put the values for the initialization
zk = [];
vk = [];
while count < itermax
    % Update x -- fmincon with 
    
    % Update z -- fmincon 
    
    % Update v -- 
    
end

%% Visualization
% TODO: what figures do we want to show?
% TODO: what are the expected results? 