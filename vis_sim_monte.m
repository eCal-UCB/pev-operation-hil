function vis_sim_monte(varargin)
% close all;
% visualizes monte carlo simulation results

if nargin == 0
    [fname, fpath] = uigetfile;
    data = load(fullfile(fpath,fname));
    sim_results = data.sim_results; 
    sim_results_base = run_sim_baseline(data); 
elseif nargin == 1
    sim_results = varargin{1}; 
    data.sim_results = sim_results;
    sim_results_base = run_sim_baseline(data); 
elseif nargin == 2
    sim_results = varargin{1}; 
    sim_results_base = varargin{2}; 
end
num_sim = length(sim_results); 
num_sim_base = length(sim_results_base);


%% Postprocessing
fprintf('[%s SIM] post-processing...\n',datetime('now')); tic;

% with controller
overstay_mean = zeros(num_sim,1);
overstay_penalty_mean = zeros(num_sim,1);
overstay_tot = zeros(num_sim,1);
profit = zeros(num_sim,1);
service_tot = zeros(num_sim,1);
for n = 1:num_sim
    overstay_mean(n) = mean(sim_results{n}.overstay_duration(sim_results{n}.overstay_duration~=0));
    overstay_penalty_mean(n) = mean(sim_results{n}.control(sim_results{n}.control(:,3)~=0,3));
    overstay_tot(n) = sum(sim_results{n}.overstay_duration);
    profit(n) = sum(sim_results{n}.profit_charging_uc+sim_results{n}.profit_charging_c+sim_results{n}.profit_overstay) - max(sim_results{n}.power) * 18.86 / 30;
    service_tot(n) = sum(sim_results{n}.num_service);
end

if num_sim_base >= 1 % if there is a baseline
    % without controller
    overstay_mean_base = zeros(num_sim_base,1);
    overstay_penalty_mean_base = zeros(num_sim_base,1);
    overstay_tot_base = zeros(num_sim_base,1);
    profit_base = zeros(num_sim_base,1);
    service_tot_base = zeros(num_sim_base,1);
    for n = 1:num_sim_base
        overstay_mean_base(n) = mean(sim_results_base{n}.overstay_duration(sim_results_base{n}.overstay_duration~=0));
        overstay_penalty_mean_base(n) = mean(sim_results_base{n}.control(sim_results_base{n}.control(:,3)~=0,3));
        overstay_tot_base(n) = sum(sim_results_base{n}.overstay_duration);
        profit_base(n) = sum(sim_results_base{n}.profit_charging_uc+sim_results_base{n}.profit_charging_c+sim_results_base{n}.profit_overstay)- max(sim_results_base{n}.power) * 18.86 / 30;
        service_tot_base(n) = sum(sim_results_base{n}.num_service);
    end
end

fprintf('[%s SIM] post-processing... DONE (%.2f sec)\n',datetime('now'),toc);


%% Visualization

% TODO: pareto chart overstay duration vs. profits with different reg params for overstaying
fprintf('[%s SIM] visualizing...\n',datetime('now')); tic;
% (1) distributions
% overstay histogram
figure('position',[1,1,640,704]); num_bins = 10; 
vals_to_vis = {'overstay_mean','profit','service_tot'};
baselines = {'overstay_mean_base', 'profit_base', 'service_tot_base'};
xlabels = {'mean overstay duration [hr]','net profit [$]','service provide [#]'};
ylabels = {'probability [0,1]'};
num_subplot = length(vals_to_vis);
for i = 1:num_subplot
    subplot(eval([num2str(num_subplot) '1' num2str(i)]));
    eval(['vals=' vals_to_vis{i} ';']); 
    eval(['baseline=' baselines{i} ';']);
    h1 = histogram(vals,num_bins,'Normalization','probability'); hold on;
    ylim([0 max(h1.BinCounts/sum(h1.BinCounts))+0.05]);
    if num_sim_base > 1
        h2 = histogram(baseline,num_bins,'Normalization','probability');
        ylim([0 max(max(h1.BinCounts/sum(h1.BinCounts)),...
            max(h2.BinCounts/sum(h2.BinCounts)))+0.05]);
    end
    s1 = stem(mean(vals),1,'b','linewidth',3, 'markersize',eps); 
    s2 = stem(mean(baseline),1,'r','linewidth',3,'markersize',eps); hold off;
    text(mean(vals),max(h1.BinCounts/sum(h1.BinCounts))+0.04, ...
        sprintf(' mean: %.1f (%+.2f%%)',mean(vals),(mean(vals)/mean(baseline)-1)*100),...
        'fontsize',15);
    grid on; xlabel(xlabels{i}); ylabel('probability $\in$ [0,1]','interpreter','latex');
    set(gca,'fontsize',15);
    if i == 1
        legend([s1 s2],{'mean w/ control','mean w/o control'},'location','northwest');
    end
end
saveas(gcf,'recent-visualization/hist.png');
saveas(gcf,'recent-visualization/hist.fig');
saveas(gcf,'recent-visualization/hist','epsc');


% (2) one day simulation results
% visualize temporal trajectories and choices
day = randi(length(sim_results));
% day = 3;
vis_sim_one_day(sim_results{day});


% (3) pareto charts
% mean overstay penalty vs. profits
figure; 
subplot(121);scatter(overstay_penalty_mean,profit,'filled'); grid on; set(gca,'fontsize',15);
xlabel('mean overstay penalty [$]'); ylabel('profit [$]');

% mean overstay penalty vs. mean overstay duration
subplot(122);scatter(overstay_penalty_mean,overstay_mean,'filled'); grid on; set(gca,'fontsize',15);
xlabel('mean overstay penalty [$]'); ylabel('mean overstay duration [hr]');
sgtitle([{'Pareto Chart'},{sprintf('(total number of simulations: %d)',length(sim_results))}],'fontsize',18);
saveas(gcf,'recent-visualization/pareto.png');
saveas(gcf,'recent-visualization/pareto.fig');
saveas(gcf,'recent-visualization/pareto','epsc');

% (4) convergence rate of BCD
ind_events = find(~isnan(sim_results{day}.choice));
max_num_iter = 0; ind_max_num_iter=0;
for e = 1:length(ind_events)
    if sim_results{day}.opts{ind_events(e)}.num_iter > max_num_iter
        max_num_iter = sim_results{day}.opts{ind_events(e)}.num_iter;
        ind_max_num_iter = ind_events(e);
    end
end
vis_sim_one_event(sim_results{day}.opts{ind_max_num_iter});

fprintf('[%s SIM] visualizing... DONE (%.2f sec)\n',datetime('now'),toc);
end