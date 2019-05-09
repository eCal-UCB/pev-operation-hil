% This script runs monte carlo simulations with various number of charging
% poles for sensitivity analysis. The number of charging poles and number
% of simulations at each monte carlo simulations are configured in 
% init_params.m. 

clear; close all;

t0=tic;
par = init_params();
fname = fullfile(pwd,'sensitivity-sim-results', ...
                sprintf('%s_senstivitiy.mat', ...
                datestr(now,'mm_dd_yy_HH_MM')));
monte_results = cell(length(par.sens_analysis.num_poles),1);
i=0;
for num_pole = par.sens_analysis.num_poles
    i = i + 1;
    save(fname);
    fprintf('\n[%s SENSITIVITY] ** start monte carlo simulation with %d poles ** \n\n',datetime('now'),num_pole);
    par.station.num_poles = num_pole;
    monte_results{i} = run_sim_monte(par); % {num_pole} -> optimal/baseline -> {num_sim}
end
save(fname);

% visualize
vis_sensitivity_analysis(monte_results);
fprintf('\n[%s SENSITIVITY] DONE sensitivity analysis, %.2f sec \n\n',datetime('now'),toc(t0));