% This script is to generate a distribution of EV charging station
% operation simulations.
%
% THIS IS A PART OF EE227C COURSE PROJECT AT UC BERKELEY.
% Contributors: Sangjae Bae, Teng Zeng, Bertrand Travacca.
% May, 2019.

close all; clear;
num_sim = 30; % simulation numbers

%% Simulation
t0 = tic;
sim_results = cell(num_sim,1);
par = set_glob_par(init_params());
if par.isFixedEventSequence
    s = 'fixed';
else
    s = 'rand';
end
fname = fullfile(pwd,'monte-sim-results',sprintf('%s_monte_sim_%d_eps_%s_seq.mat',datestr(now,'mm_dd_yy_HH_MM'),num_sim,s));
for n = 1:num_sim
    save(fname);
    fprintf('======================== %d/%d =======================\n',...
        n,num_sim);
    t1= tic;
    sim_results{n} = run_sim_one_day(par);
    fprintf('\n[%s SIM] one day operation DONE (%.2f sec)\n\n\n',datetime('now'),toc(t1));
end
tot_time = toc(t0);
fprintf('[%s SIM] total computation time: %.2f sec\n',datetime('now'), tot_time);
save(fname);

%% Visualization
vis_sim_monte(sim_results);
