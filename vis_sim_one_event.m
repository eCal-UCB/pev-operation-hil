function vis_sim_one_event(opt)
xticklabel_ = {'Flex, z_{c}';'ASAP, z_{uc}';'Overstay, y'};

% simulation result
figure('position',[1,1,640,704]);
subplot(411)
plot(opt.prb.user.time:opt.par.Ts:opt.prb.user.time+opt.prb.user.duration-opt.par.Ts, opt.x(opt.prb.N_flex+2:end),'linewidth',1.5); grid on;
xlabel('opt.parking durations (hour)'); ylabel('kW'); set(gca,'fontsize',15);
subplot(412)
plot(opt.prb.user.time:opt.par.Ts:opt.prb.user.time+opt.prb.user.duration, opt.x(1:opt.prb.N_flex+1),'linewidth',1.5); grid on;
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
% plot(opt.J(opt.J~=0),'linewidth',1.5);grid on;hold on;
% plot(1:length(opt.J(opt.J~=0)),opt.J(length(opt.J(opt.J~=0)))*ones(size(opt.J(opt.J~=0))),'--r','linewidth',1.5); hold off;
plot(log(opt.J(opt.J~=0)/min(opt.J(opt.J~=0))),'linewidth',1.5);grid on;hold on;
plot(1:length(opt.J(opt.J~=0)),zeros(size(opt.J(opt.J~=0))),'--r','linewidth',1.5); hold off;
xlim([1,length(opt.J(opt.J~=0))]);
xlabel('iteration'); 
ylabel({'$log(\frac{cost}{cost^\ast})$'},'interpreter','latex'); 
title([{'BCD Convergence'},sprintf('(total iter: %d)',opt.num_iter)]);
set(gca,'fontsize',15);
saveas(gcf,'recent-visualization/convergence.png');
saveas(gcf,'recent-visualization/convergence.fig');
saveas(gcf,'recent-visualization/convergence','epsc');
end