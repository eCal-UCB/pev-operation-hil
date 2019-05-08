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
J = @(z,x,v) dot([(sum((x(prb.N_flex+2:end).*(prb.TOU(1:prb.N_flex) - z(1)))+par.lambda.x.*x(prb.N_flex+2:end))+par.lambda.z_c*z(1)^2) + par.lambda.h_c * 1/z(3);
            sum((prb.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2)))+par.lambda.z_uc*z(2)^2) + par.lambda.h_uc * 1/z(3);
            1/3*sum(prb.station.pow_max*(prb.TOU(1:prb.N_asap) - 0))],v);

        
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
    if mod(count,1) == 0
        fprintf('[ OPT] iter: %d, Jk: %.3f\n',count,Jk(count));
    end
    
    % update init variables
    prb.z0 = zk; prb.x0 = xk; prb.v0 = vk; set_glob_prb(prb);
    
    % update control variables
    zk = argmin_z([],xk,vk);
    xk = argmin_x(zk,[],vk);
    vk = argmin_v(zk,xk,[]);
    
    % compute residual
    improve = Jk(count)-J(zk,xk,vk);
end
opt.z = zk;
opt.tariff.flex = zk(1);
opt.tariff.asap = zk(2);
opt.tariff.overstay = zk(3);
opt.x = xk;
opt.flex.SOCs = xk(1:prb.N_flex+1);
opt.flex.powers = xk(prb.N_flex+2:end);
opt.asap.powers = prb.station.pow_max;
opt.v = vk;
opt.prob.flex = vk(1);
opt.prob.asap = vk(2);
opt.prob.leave = vk(3);
opt.J = Jk(1:count);
opt.num_iter = count;
opt.prb = prb;
opt.par = par;
opt.time.start = prb.user.time;
opt.time.end_flex = prb.user.time + prb.user.duration;
opt.time.end_asap = prb.user.time + prb.N_asap*par.Ts;
opt.isOverStay = false;
fprintf('[ OPT] DONE (%.2f sec) sum(vk) = %.2f, iterations = %d\n',toc,sum(vk),count);


%% Visualization
% TODO: what figures do we want to show?
% TODO: what are the expected results? 
disp('[ VIS] visualizing...');tic;
xticklabel_ = {'Flex, z_{c}';'ASAP, z_{uc}';'Overstay, y'};

% simulation result
figure('position',[1,1,640,704]);
subplot(411)
plot(opt.prb.user.time:opt.par.Ts:opt.prb.user.time+opt.prb.user.duration-opt.par.Ts, opt.x(opt.prb.N_flex+2:end),'linewidth',1.5); grid on;
xlim([prb.user.time prb.user.time+prb.user.duration])
xlabel('opt.parking durations (hour)'); ylabel('kW'); set(gca,'fontsize',15);
subplot(412)
plot(opt.prb.user.time:opt.par.Ts:opt.prb.user.time+opt.prb.user.duration, opt.x(1:opt.prb.N_flex+1),'linewidth',1.5); grid on;
xlim([prb.user.time prb.user.time+prb.user.duration])
xlabel('opt.parking durations (hour)'); ylabel('SOC [0,1]'); set(gca,'fontsize',15);
subplot(413)
bar(opt.z(1:3)); grid on; ylabel('$/kW');
set(gca,'fontsize',15,'xticklabel',xticklabel_);
subplot(414)
bar(opt.v); grid on; ylabel('Probability')
set(gca,'fontsize',15,'xticklabel',xticklabel_);
saveas(gcf,'recent-visualization/one_event.png');
saveas(gcf,'recent-visualization/one_event.fig');
saveas(gcf,'recent-visualization/one_event','epsc');


% convergence analysis
figure;
plot(opt.J(opt.J~=0),'linewidth',1.5);grid on;hold on;
plot(1:length(opt.J(opt.J~=0)),opt.J(length(opt.J(opt.J~=0)))*ones(size(opt.J(opt.J~=0))),'--r','linewidth',1.5); hold off;
xlim([1,length(opt.J(opt.J~=0))]);
xlabel('Iteration'); ylabel('Cost'); title([{'BCD Convergence'},sprintf('(total iter: %d)',opt.num_iter)]);
set(gca,'fontsize',15);
saveas(gcf,'recent-visualization/convergence.png');
saveas(gcf,'recent-visualization/convergence.fig');
saveas(gcf,'recent-visualization/convergence','epsc');

fprintf('[ VIS] DONE (%.2f sec)\n',toc);