% This script is to debug the station wide optimziation problem by
% comparing with the single charger optimization. 
%
% Example1:
% - For one user, station-opt must return the same optimal solution
% with single-opt
%
% Example 2:
% - For two users, station-opt must return optimal solution that has
% higher expected revenue than single-opt


%% Generate event
par = set_glob_par(init_params());
events = gen_events_one_day(par);
event = events.inp{1}; % fix the event with the first one
set_glob_prb(init_prb(event));

%% Example 1
% Run single-opt
opt_single = run_opt();

% Run station-opt
station = init_station(par);
[station, opt_station] = run_opt_station(station,  round(event.time/par.Ts)*par.Ts);

% compute the cost
opt_single.cost = sum((opt_single.prb.TOU(1:opt_single.prb.N_flex)-opt_single.z(1)).*opt_single.x(1:opt_single.prb.N_flex)) + station('cost_dc')*(opt_single.peak_pow);
opt_station.cost = sum((opt_station.prb.TOU(1:opt_station.prb.N_flex)-opt_station.z(1)).*opt_station.x(1:opt_station.prb.N_flex)) + station('cost_dc')*(opt_station.peak_pow);

%%
fprintf('single: %.1f, station: %.1f\n',opt_single.cost, opt_station.cost)