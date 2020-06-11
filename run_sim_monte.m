function varargout = run_sim_monte(varargin)
% This script is to generate a distribution of EV charging station
% operation simulations.
%
% THIS IS A PART OF EE227C COURSE PROJECT AT UC BERKELEY.
% Contributors: Sangjae Bae, Teng Zeng, Bertrand Travacca.
% May, 2019.

% close all; clear;

%% Simulation
t0 = tic;

if nargin == 0
    par = set_glob_par(init_params());
elseif nargin == 1
    par = set_glob_par(varargin{1});
end

if par.sim.isFixedEventSequence
    s = 'fixed';
else
    s = 'rand';
end

num_sim = par.monte.num_sims; % simulation numbers
sim_results = cell(num_sim,1);
sim_results_base = cell(num_sim,1);
fname = fullfile(pwd,'monte-sim-results', ...
                sprintf('%s_monte_eps%d_%s_seq_poles%d.mat', ...
                datestr(now,'mm_dd_yy_HH_MM'),num_sim,s,par.station.num_poles));
for n = 1:num_sim
    save(fname);
    fprintf('======================== %d/%d =======================\n',...
        n,num_sim);
    t1= tic;
%     sim_results{n} = run_sim_one_day(par);
    sim_results{n} = run_sim_one_day_v2(par);
    sim_results_base{n} = run_sim_one_day_baseline(sim_results{n});
    fprintf('\n[%s SIM] one day operation DONE (%.2f sec)\n\n\n',datetime('now'),toc(t1));
end
tot_time = toc(t0);
fprintf('[%s SIM] DONE monte carlo simulation, %.2f sec\n',datetime('now'), tot_time);
save(fname);

if nargout == 1 
    varargout = {};
    arg = {};
    arg.optimal = sim_results;
    arg.baseline = sim_results_base;
    varargout{1} = arg;
end

%% Visualization
if nargin == 0
    vis_sim_monte(sim_results,sim_results_base);
    vis_sim_monte_vk(sim_results,sim_results_base);
end

end
