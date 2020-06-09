function varargout = run_sim_one_day_v2(varargin)
% clear;clc;close all;
% This script is to simulate EV charging station operations where the
% charging tariff is determined real-time with taking account into EV
% drivers' behaviors. The overall objective of the tariff control is to
% minimize the operation cost of the system operater. 
%
% There are three choices that each EV driver can make at arrival:
% (i) charging with flexibility, (ii) charging as soon as possible, and
% (iii) leaving without charging.
%
% At each arrival of EV driver, user specific parameters, e.g., battery
% capacity, desired parking durations, initial SOC, and needed SOC level
% (for next mobility demand) are randomly sampled from an empirical
% probability distribution function that is generated with TELD dataset 
% [X].
%
% last modified, March 2020 - Teng
% v2 - entire station optimization
% v1 - single charger optimization
%
% Contributors: Sangjae Bae, Teng Zeng, Bertrand Travacca.

% clear

%% Initialization
fprintf('[%s INIT] initializing...\n',datetime('now'));
if nargin == 0
    par = set_glob_par(init_params());
    events = gen_events_one_day(par);
elseif nargin == 1
    par = varargin{1};
    events = gen_events_one_day(par);
elseif nargin == 2
    par = varargin{1};
    events = varargin{2};
else
    error(sprintf('[%s ERROR] too many input arguments',datetime('now')));
end

fprintf('[%s INIT] DONE\n',datetime('now'));

%% Simulation
t = par.sim.starttime:par.Ts:par.sim.endtime; i_k = 0; i_event = 0;
sim = init_sim(t); % simulation result
station = containers.Map; % station monitor
station('num_occupied_pole') = 0; 
station('FLEX_list') = [];
station('ASAP_list') = [];
station('num_empty_pole') = par.station.num_poles;
station('D_init') = 0;
station('pow_cap') = 6*6.6 ; % this value is arbitrary for now
station('cost_dc') = 20; % this value is arbitrary for now
sim.events = events;

for k = par.sim.starttime:par.Ts:par.sim.endtime
    i_k = i_k + 1;
    % check visit
    if i_event <= length(events.time)
        if any(round(events.time/par.Ts)*par.Ts == k)
            inds_events_k = find(round(events.time/par.Ts)*par.Ts == k);
            for j = 1:length(inds_events_k)
                i_event = i_event + 1; % number of investigated events
                if events.inp{inds_events_k(j)}.duration <= par.sim.endtime - k ...
                        && station('num_empty_pole') > 0
                   sim.tot_decision = sim.tot_decision + 1; % number of decisions
                   sim.events.triggered(i_event) = true; % this event is triggered
                   
                   set_glob_prb(init_prb(events.inp{inds_events_k(j)}));

                   % find optimal tariff
                   if isempty(station('FLEX_list')) && isempty(station('ASAP_list'))
                       opt = run_opt();
                   else
                       [station, opt] = run_opt_station(station, k);
                   end
                   sim.opts{i_event} = opt;
                   station('D_init') = opt.peak_pow; % update demand charge
                   
                   % driver makes choice
                   if par.sim.isFixedSeed
                       rng(1);
                   end
                   rc = rand;
                   if rc <= opt.prob.flex
                       opt.choice = 0; % charging with flexibility
                       opt.time.end = opt.time.end_flex;
                       opt.powers = opt.flex.powers; % TODO: need to adjust the power profile for existing users with right time index
                       opt.price = opt.tariff.flex;
                       if isempty(station('FLEX_list')) % record FLEX user
                           station('FLEX_list') = {['EV' num2str(sim.tot_decision)]};
                       else
                           station('FLEX_list') = [station('FLEX_list') {['EV' num2str(sim.tot_decision)]}];
                       end
                   elseif rc <= opt.prob.flex + opt.prob.asap
                       opt.choice = 1; % charging as soon as possible
                       opt.time.end = opt.time.end_asap;
                       opt.powers = opt.asap.powers;
                       opt.price = opt.tariff.asap;
                       if isempty(station('ASAP_list')) % record FLEX user
                           station('ASAP_list') = {['EV' num2str(sim.tot_decision)]};
                       else
                           station('ASAP_list') = [station('ASAP_list') {['EV' num2str(sim.tot_decision)]}];
                       end
                   else
                       opt.choice = 2; % leaving without charging
                   end
                   sim.choice_probs(i_event,:) = opt.v;
                   sim.choice(i_event) = opt.choice;
                   sim.control(i_event,:) = opt.z(1:3);
                   fprintf('[%s EVENT] time = %.2f, CHOICE = %s\n',datetime('now'),k,par.dcm.choices{opt.choice+1});

                   % if the driver chooses to charge EV
                   if opt.choice <= 1
                       [opt.time.leave, duration] = get_rand_os_duration(opt);
                       try
                           sim.overstay_duration(i_k) = sim.overstay_duration(i_k) + duration;
                           if duration == 32
                               disp('--- duration is boomed ----')
                           end
                       catch
                           a = 1;
                       end
                       sim.num_service(i_k) = sim.num_service(i_k) + 1;
                       station('num_occupied_pole') = station('num_occupied_pole') + 1;
                       station('num_empty_pole') = station('num_empty_pole') - 1;
                       opt.power_traj_actual = [];
                       station(['EV' num2str(sim.tot_decision)]) = opt;
                   end 
                else
                    if station('num_empty_pole') == 0
                        fprintf('[%s EVENT] SKIPPED (event %d) due to full occupancy\n',datetime('now'),i_event);
                    else
                        fprintf('[%s EVENT] SKIPPED (event %d) due to violating operationg hours\n',datetime('now'),i_event);
                    end
                end
            end
        end
    end
    
    % update agg
    keys = station.keys();
    if ~isempty(keys)
        for ev = keys
            if contains(ev{1},'EV')
                if  k < station(ev{1}).time.end % is charging duration
                    TOU = interp1(0:0.25:24-0.25,par.TOU,k,'nearest');
                    % add actual power_ record to user
                    opt = station(ev{1});
                    power = station(ev{1}).powers(1);
                    opt.power_traj_actual = [opt.power_traj_actual power];
                    station(ev{1}) = opt;
%                     if length(station(ev{1}).powers) > 1
%                         power = interp1(linspace(station(ev{1}).time.start, ...
%                                             station(ev{1}).time.end,...
%                                             length(station(ev{1}).powers)), ...
%                                             station(ev{1}).powers, k);
%                     elseif length(station(ev{1}).powers) == 1
%                         power = station(ev{1}).powers;
%                     end
                    if power == opt.prb.station.pow_max % hyperthetically when power is max power it's the uncontrol charging
                        sim.profit_charging_uc(i_k) = sim.profit_charging_uc(i_k) + par.Ts * power * (station(ev{1}).price - TOU);
                    else % flexible charging
                        sim.profit_charging_c(i_k) = sim.profit_charging_c(i_k) + par.Ts * power * (station(ev{1}).price - TOU);
                    end
                    sim.power(i_k) = sim.power(i_k) + power;
                    
%                     sim.profit_charging(i_k) = sim.profit_charging(i_k) + par.Ts * power * (station(ev{1}).price - TOU);
                    sim.occ.charging(i_k) = sim.occ.charging(i_k) + 1;
                else % is overstaying
                    if k <= station(ev{1}).time.leave 
                        sim.profit_overstay(i_k) = sim.profit_overstay(i_k) + par.Ts * station(ev{1}).tariff.overstay;
                        sim.occ.overstay(i_k) = sim.occ.overstay(i_k) + 1;
%                         sim.overstay_duration(i_k) = sim.overstay_duration(i_k) + par.Ts;
                    else
                        station.remove(ev{1});
                            
                        station('num_occupied_pole') = station('num_occupied_pole') - 1;
                        station('num_empty_pole') = station('num_empty_pole') + 1;
                    end 
                    
                    % remove ev from flex and asap list, it needs no
                    % charging
                    flex_list = station('FLEX_list');
                    asap_list = station('ASAP_list');
                    if any(strcmp(flex_list,ev{1})) % remove from flex list
                        rm_idx = find(contains(flex_list,ev{1}));
                        station('FLEX_list') = [flex_list(1:rm_idx-1) flex_list(rm_idx+1:end)];
                    elseif any(strcmp(asap_list,ev{1})) % remove from asap list
                        rm_idx = find(contains(asap_list,ev{1}));
                        station('ASAP_list') = [asap_list(1:rm_idx-1) asap_list(rm_idx+1:end)];
                    end
                end
            elseif contains(ev{1},'occ')
                sim.occ.total(i_k) = station('num_occupied_pole');
                sim.occ.empty(i_k) = station('num_empty_pole');
            end
        end
    end
end

sim.par = par;

varargout = {};
varargout{1} = sim;

%% Visualization
if nargout == 0
    options = vis_sim_one_day(); % options: display, temporals, choices
    options.temporals = true;
    vis_sim_one_day(sim,options);
end

end