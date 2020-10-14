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
sim_results_station = cell(num_sim,1);
sim_results_base = cell(num_sim,1);
fname = fullfile('../monte-sim-results', ...
                sprintf('%s_monte_eps%d_%s_seq_poles%d.mat', ...
                datestr(now,'mm_dd_yy_HH_MM'),num_sim,s,par.station.num_poles));
n = 1;

if par.sim.isFixedSeed
    seed_list = linspace(1, num_sim, num_sim);
end
while n <= num_sim
            % for n = 1:num_sim
    save(fname);
    fprintf('======================== %d/%d =======================\n',...
        n,num_sim);
    t1= tic;
    try
        if par.sim.isFixedSeed
            seed_val = seed_list(n);
            events = gen_events_one_day(par, seed_val);
        else
            events = gen_events_one_day(par);
        end
        disp('Optimize single charger...');
        sim_results{n} = run_sim_one_day(par,events);
        disp('Optimize station...');
        sim_results_station{n} = run_sim_one_day_station(par,events);
        disp('Run baseline...');
        sim_results_base{n} = run_sim_one_day_baseline(sim_results{n});
        fprintf('\n[%s SIM] one day operation DONE (%.2f sec)\n\n\n',datetime('now'),toc(t1));
        n = n + 1;
    catch
        warning(sprintf('Constraint violated. Regenerate count %d',n)); %#ok<SPWRN>
    end
end
tot_time = toc(t0);
fprintf('[%s SIM] DONE monte carlo simulation, %.2f sec\n',datetime('now'), tot_time);
save(fname);

if nargout == 1 
    varargout = {};
    arg = {};
    arg.optimal = sim_results;
    arg.optimal_station = sim_results_station;
    arg.baseline = sim_results_base;
    varargout{1} = arg;
end

%% Visualization
if nargin == 0
    vis_sim_monte(sim_results,sim_results_base);
    vis_sim_monte(sim_results_station,sim_results_base);
    vis_sim_monte(sim_results_station,sim_results);
    vis_sim_monte_vk(sim_results,sim_results_base);
end

end
