% This script is to generate a distribution of EV charging station
% operation simulations.
%
% THIS IS A PART OF EE227C COURSE PROJECT AT UC BERKELEY.
% Contributors: Sangjae Bae, Teng Zeng, Bertrand Travacca.
% May, 2019.

close all; clear;
num_sim = 50; % simulation numbers

%% Simulation
t0 = tic;
sim_results = cell(num_sim,1);
fname = fullfile(pwd,'monte-sim-results',sprintf('%s monte-sim-%d-eps.mat',datetime('now'),num_sim));
for n = 1:num_sim
    save(fname);
    fprintf('======================== %d/%d =======================\n',...
        n,num_sim);
    t1= tic;
    sim_results{n} = run_sim_one_day();
    fprintf('\n[%s SIM] one day operation DONE (%.2f sec)\n\n\n',datetime('now'),toc(t1));
end
tot_time = toc(t0);
fprintf('[%s SIM] total computation time: %.2f sec\n',datetime('now'), tot_time);
save(fname);

%% Visualization
vis_sim_monte(sim_results);
