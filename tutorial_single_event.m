% This script is to simulate an instance of an EV charging tariff 
% controller that controls charging tariff of two charging options:
% (i) Charging with flexibitliy
% (ii) Charging as soon as possible
% 
% EE227C project, May 2019.

clear; tic;

%% Initialize parameters
disp('[ INIT] initializing...');
par = set_glob_par(init_params());
prb = set_glob_prb(init_prb());
disp('[ INIT] DONE');


%% Functions
J = @(z,x,v) dot([sum((x(prb.N_flex+2:end).*(prb.TOU(1:prb.N_flex) - z(1))).^2) + par.lambda.h_c * 1/z(3); % h_c
            sum((par.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2))).^2) + par.lambda.h_uc * 1/z(3); % h_uc
            sum((par.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2))).^2)],v); % h_l

        
%% Run algorithm -- block coordinate descent
itermax = 1e4;
count = 0; improve = inf;
zk = ones(4,1);                         % [z_c, z_uc, y, 1];
xk = ones(2*prb.N_flex+1,1);            % [soc0, ..., socN, u0, ..., uNm1];
vk = 1/3*ones(3,1);                     % [sm_c, sm_uc, sm_y];
Jk = zeros(itermax,1);
while count < itermax && improve >= 0 && abs(improve) >= par.opt.eps
    count = count + 1;
    Jk(count) = J(zk,xk,vk);    
    
    % update init variables
    prb.z0 = zk; prb.x0 = xk; prb.v0 = vk; set_glob_prb(prb);
    
    % update control variables
    zk = argmin_z([],xk,vk);
    xk = argmin_x(zk,[],vk);
    vk = argmin_v(zk,xk,[]);
    
    % compute residual
    improve = Jk(count)-J(zk,xk,vk);
    
    if mod(count,1) == 0
        fprintf('[ OPT] iter: %d, improve: %.3f\n',count,improve);
    end
end

fprintf('[ OPT] DONE (%.2f sec) sum(vk) = %.2f, iterations = %d\n',toc,sum(vk),count);


%% Visualization
% TODO: what figures do we want to show?
% TODO: what are the expected results? 
disp('[ VIS] visualizing...');tic;
xticklabel_ = {'Flex, z_{c}';'ASAP, z_{uc}';'Overstay, y'};

% simulation result
figure(1)
subplot(411)
plot(prb.user.time:par.Ts:prb.user.time+prb.user.duration-par.Ts, xk(prb.N_flex+2:end),'linewidth',1.5); grid on;
xlabel('Parking durations (hour)'); ylabel('kW'); set(gca,'fontsize',15);
subplot(412)
plot(prb.user.time:par.Ts:prb.user.time+prb.user.duration, xk(1:prb.N_flex+1),'linewidth',1.5); grid on;
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
xlim([1,length(Jk(Jk~=0))]);
xlabel('Iteration'); ylabel('Cost'); title('BCD Convergence');
set(gca,'fontsize',15);

fprintf('[ VIS] DONE (%.2f sec)\n',toc);