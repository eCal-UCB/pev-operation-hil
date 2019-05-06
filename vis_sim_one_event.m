function vis_sim_one_event(opt)
xticklabel_ = {'Flex, z_{c}';'ASAP, z_{uc}';'Overstay, y'};

% simulation result
figure;
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

% convergence analysis
figure;
plot(opt.J(opt.J~=0),'linewidth',1.5);grid on;hold on;
plot(1:length(opt.J(opt.J~=0)),opt.J(length(opt.J(opt.J~=0)))*ones(size(opt.J(opt.J~=0))),'--r','linewidth',1.5); hold off;
xlim([1,length(opt.J(opt.J~=0))]);
xlabel('Iteration'); ylabel('Cost'); title([{'BCD Convergence'},sprintf('(total iter: %d)',opt.num_iter)]);
set(gca,'fontsize',15);
end