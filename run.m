% This script is to run a simulation of an EV charging tariff controller
% that controls charging tariff of two charging options:
% (i) Charging with flexibitliy
% (ii) Charging as soon as possible
%
% EE227C project, May 2019.

clear; tic;

%% Initialize parameters
set_glob_par(init_params());
par = get_glob_par();


%% Helper functions
J = @(z,x,v) dot([sum((x(par.N_flex+2:end).*(par.TOU(1:par.N_flex*par.Ts) - z(1))).^2) + par.lambda.h_c * 1/z(3); % h_c
            sum((par.station.pow_max*(par.TOU(1:par.N_asap*par.Ts) - z(2))).^2) + par.lambda.h_uc * 1/z(3); % h_uc
            sum((par.station.pow_max*(par.TOU(1:par.N_asap*par.Ts) - z(2))).^2)],v); % h_l

        
%% Run algorithm -- block coordinate descent
itermax = 1e4;
count = 0; residual = 100;
zk = ones(4,1);                         % [z_c, z_uc, y, 1];
xk = ones(2*par.N_flex+1,1);            % [soc0, ..., socN, u0, ..., uNm1];
vk = 1/3*ones(3,1);                     % [sm_c, sm_uc, sm_y];
while count < itermax && residual > 1
    J0=J(zk,xk,vk);
    par.z0 = zk; par.x0 = xk; par.v0 = vk; set_glob_par(par);
    
    % update
    zk = argmin_z([],xk,vk);
    xk = argmin_x(zk,[],vk);
    vk = argmin_v(zk,xk,[]);
    
    residual = abs(J(zk,xk,vk)-J0);
    count = count + 1;
    if mod(count,10) == 0
        fprintf('iter: %d, J: %.2f\n',count,J(zk,xk,vk));
    end
end

fprintf('[ DONE] (%.2f sec) sum(vk) = %.2f, iterations = %d\n',toc,sum(vk),count);


%% Visualization
% TODO: what figures do we want to show?
% TODO: what are the expected results? 

% simulation result
figure(1)
subplot(411)
plot(xk(par.N_flex+2:end),'linewidth',1.5); grid on;
xlabel('Parking durations (hour)'); ylabel('kW'); set(gca,'fontsize',15);
subplot(412)
plot(xk(1:par.N_flex+1),'linewidth',1.5); grid on;
xlabel('Parking durations (hour)'); ylabel('SOC [0,1]'); set(gca,'fontsize',15);
subplot(413)
bar(zk(1:3)); grid on; ylabel('$/kW');
set(gca,'fontsize',15,'xticklabel',{'Flex charging, z_{c}';'ASAP charging z_{uc}';'Overstay penalty, y'});
subplot(414)
bar(vk); grid on; ylabel('Probability')
set(gca,'fontsize',15,'xticklabel',{'Flex charging, z_{c}';'ASAP charging z_{uc}';'Overstay penalty, y'});

% convergence analysis
% figure(2)