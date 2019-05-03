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
J = @(z,x,v) dot([sum((x(par.N_flex+2:end).*(par.TOU(1:par.N_flex) - z(1))).^2) + par.lambda.h_c * 1/z(3)^2; % h_c
            sum((par.station.pow_max*(par.TOU(1:par.N_asap) - z(2))).^2) + par.lambda.h_uc * 1/z(3)^2; % h_uc
            sum((par.station.pow_max*(par.TOU(1:par.N_asap) - z(2))).^2)],v); % h_l

        
%% Run algorithm -- block coordinate descent
itermax = 1e4;
count = 1; residual = 100;
zk = ones(4,1);                         % [z_c, z_uc, y, 1];
xk = ones(2*par.N_flex+1,1);            % [soc0, ..., socN, u0, ..., uNm1];
vk = 1/3*ones(3,1);                     % [sm_c, sm_uc, sm_y];
Jk = zeros(itermax,1);
while count < itermax && residual > par.opt.eps
    Jk(count) = J(zk,xk,vk);
    if mod(count,10) == 0
        fprintf('[ SOLVING] iter: %d, cost: %.2f\n',count,Jk(count));
    end
    
    % update init variables
%     J0=J(zk,xk,vk);
    par.z0 = zk; par.x0 = xk; par.v0 = vk; set_glob_par(par);
    
    % update control variables
    zk = argmin_z([],xk,vk);
    xk = argmin_x(zk,[],vk);
    vk = argmin_v(zk,xk,[]);
    
    % update condition variables
    residual = abs(J(zk,xk,vk)-Jk(count));
    count = count + 1;
end

fprintf('[ DONE] (%.2f sec) sum(vk) = %.2f, iterations = %d\n',toc,sum(vk),count);


%% Visualization
% TODO: what figures do we want to show?
% TODO: what are the expected results? 
xticklabel_ = {'Flex, z_{c}';'ASAP, z_{uc}';'Overstay, y'};

% simulation result
figure(1)
subplot(411)
plot(0:par.Ts:par.user.duration-par.Ts, xk(par.N_flex+2:end),'linewidth',1.5); grid on;
xlabel('Parking durations (hour)'); ylabel('kW'); set(gca,'fontsize',15);
subplot(412)
plot(0:par.Ts:par.user.duration, xk(1:par.N_flex+1),'linewidth',1.5); grid on;
xlabel('Parking durations (hour)'); ylabel('SOC [0,1]'); set(gca,'fontsize',15);
subplot(413)
bar(zk(1:3)); grid on; ylabel('$/kW');
set(gca,'fontsize',15,'xticklabel',xticklabel_);
subplot(414)
bar(vk); grid on; ylabel('Probability')
set(gca,'fontsize',15,'xticklabel',xticklabel_);

% convergence analysis
figure(2)
plot(Jk(Jk~=0),'linewidth',1.5);grid on;hold on;
plot(1:length(Jk(Jk~=0)),Jk(length(Jk(Jk~=0)))*ones(size(Jk(Jk~=0))),'--r','linewidth',1.5); hold off;
xlabel('Iteration'); ylabel('Cost'); title('BCD Convergence');
set(gca,'fontsize',15);
