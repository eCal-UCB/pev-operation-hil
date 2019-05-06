% This script is to generate a distribution of operation simulation results
% for optimal EV charging station control. 

close all; clear;
num_sim = 10;


%% Simulation
t0 = tic;
sim_results = cell(num_sim,1);
for n = 1:num_sim
    fprintf('======================== %d/%d =======================\n',...
        n,num_sim);
    t1= tic;
    sim_results{n} = run_sim_one_day();
    fprintf('\n[%s SIM] one day operation DONE (%.2f sec)\n\n\n',datetime('now'),toc(t1));
end
tot_time = toc(t0);
fprintf('[%s SIM] total computation time: %.2f sec\n',datetime('now'), tot_time);
save(fullfile(pwd,'monte-sim-results',sprintf('%s monte-sim.mat',datetime('now'))));


%% Postprocessing
fprintf('[%s SIM] post-processing...\n',datetime('now')); tic;
overstay_mean = zeros(num_sim,1);
overstay_tot = zeros(num_sim,1);
profit = zeros(num_sim,1);
overstay_penalty_mean = zeros(num_sim,1);
for n = 1:num_sim
    overstay_mean(n) = mean(sim_results{n}.overstay_duration(sim_results{n}.overstay_duration~=0));
    overstay_tot(n) = sum(sim_results{n}.overstay_duration);
    profit(n) = sum(sim_results{n}.profit);
    overstay_penalty_mean(n) = mean(sim_results{n}.control(sim_results{n}.control(:,3)~=0,3));
end
fprintf('[%s SIM] post-processing... DONE (%.2f)\n',datetime('now'),toc);


%% Visualization
% TODO: pareto chart overstay duration vs. profits with different reg params for overstaying

% (1) distributions
% overstay histogram
figure(1); num_bins = 10;
subplot(211);
h1 = histogram(overstay_mean,num_bins,'Normalization','probability'); hold on;
s1 = stem(mean(overstay_mean),1,'b','linewidth',3, 'markersize',eps); hold off;
ylim([0 max(h1.BinCounts/sum(h1.BinCounts))+0.05]);
text(mean(overstay_mean),max(h1.BinCounts/sum(h1.BinCounts))+0.04,sprintf(' mean: %.2f hours',mean(overstay_mean)),'fontsize',15);
grid on; xlabel('average overstay duration'); ylabel('probability [0,1]');
set(gca,'fontsize',15); legend(s1,{'mean'},'location','northwest');

% profit histogram
subplot(212);
h2 = histogram(profit,num_bins,'Normalization','probability'); hold on;
stem(mean(profit),1,'b','linewidth',3); hold off;
ylim([0 max(h2.BinCounts/sum(h2.BinCounts))+0.05]);
text(mean(profit),max(h2.BinCounts/sum(h2.BinCounts))+0.04,sprintf(' mean: $ %.2f',mean(profit)),'fontsize',15);
grid on; xlabel('total profit'); ylabel('probability [0,1]');
set(gca, 'fontsize', 15);


% (2) one day simulation results
% visualize temporal trajectories and choices
day = randi(length(sim_results));
vis_sim_one_day(sim_results{day});


% (3) pareto charts
% mean overstay penalty vs. profits
figure; 
subplot(121);scatter(overstay_penalty_mean,profit,'filled'); grid on; set(gca,'fontsize',15);
xlabel('mean overstay penalty ($)'); ylabel('profit ($)');

% mean overstay penalty vs. mean overstay duration
subplot(122);scatter(overstay_penalty_mean,overstay_mean,'filled'); grid on; set(gca,'fontsize',15);
xlabel('mean overstay penalty ($)'); ylabel('mean overstay duration (hour)');
sgtitle([{'Pareto Chart'},{sprintf('(total simulation num: %d)',length(sim_results))}],'fontsize',18);


% (4) convergence rate of BCD
ind_events = find(~isnan(sim_results{day}.choice));
max_num_iter = 0; ind_max_num_iter=0;
for e = ind_events
    if sim_results{day}.opts{e}.num_iter > max_num_iter
        max_num_iter = sim_results{day}.opts{e}.num_iter;
        ind_max_num_iter = e;
    end
end
vis_sim_one_event(sim_results{day}.opts{ind_max_num_iter});