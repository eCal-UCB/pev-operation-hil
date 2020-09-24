% This script runs power management optimization only.
function varargout = run_power_management(varargin)

%% Initialization
if nargin == 0
    par = set_glob_par(init_params());
    if par.sim.isFixedSeed
        seed_val = 3;
        events = gen_events_one_day(par, seed_val);
    else
        events = gen_events_one_day(par);
    end
elseif nargin == 2
    par = varargin{1};
    events = varargin{2};
end


%% Simulation -- single charger
t = par.sim.starttime:par.Ts:par.sim.endtime; i_k = 0; i_event = 0;
sim_single = init_sim(t); % simulation result
single = containers.Map; % station monitor
single('num_occupied_pole') = 0; 
single('num_empty_pole') = par.station.num_poles;
sim_single.events = events;

for k = par.sim.starttime:par.Ts:par.sim.endtime
    i_k = i_k + 1;
    % check visit
    if i_event <= length(events.time)
        if any(round(events.time/par.Ts)*par.Ts == k)
            inds_events_k = find(round(events.time/par.Ts)*par.Ts == k);
            for j = 1:length(inds_events_k)
                i_event = i_event + 1; % number of investigated events
                if events.inp{inds_events_k(j)}.duration <= par.sim.endtime - k ...
                        && single('num_empty_pole') > 0
                   sim_single.tot_decision = sim_single.tot_decision + 1; % number of decisions
                   sim_single.events.triggered(i_event) = true; % this event is triggered
                    
                   set_glob_prb(init_prb(events.inp{inds_events_k(j)}));

                   % find optimal tariff
                   opt = run_opt();
                   sim_single.opts{i_event} = opt;

                   % driver makes choice
                   rc = 0; % XXX always choose flexible.
                   if rc <= opt.prob.flex
                       opt.choice = 0; % charging with flexibility
                       opt.time.end = opt.time.end_flex;
                       opt.powers = opt.flex.powers;
                       opt.price = opt.tariff.flex;
                   elseif rc <= opt.prob.flex + opt.prob.asap
                       opt.choice = 1; % charging as soon as possible
                       opt.time.end = opt.time.end_asap;
                       opt.powers = opt.asap.powers;
                       opt.price = opt.tariff.asap;
                   else
                       opt.choice = 2; % leaving without charging
                   end
                   sim_single.choice_probs(i_event,:) = opt.v;
                   sim_single.choice(i_event) = opt.choice;
                   sim_single.control(i_event,:) = opt.z(1:3);
                   if par.VIS_DETAIL
                    fprintf('[%s EVENT] time = %.2f, CHOICE = %s\n',datetime('now'),k,par.dcm.choices{opt.choice+1});
                   end

                   % if the driver chooses to charge EV
                   if opt.choice <= 1
                       [opt.time.leave, duration] = get_rand_os_duration(opt);
                       sim_single.overstay_duration(i_k) = sim_single.overstay_duration(i_k) + duration;
                       sim_single.num_service(i_k) = sim_single.num_service(i_k) + 1;
                       single('num_occupied_pole') = single('num_occupied_pole') + 1;
                       single('num_empty_pole') = single('num_empty_pole') - 1;
                       single(['EV' num2str(sim_single.tot_decision)]) = opt;
                   end 
                else
                    if par.VIS_DETAIL
                        if single('num_empty_pole') == 0
                            fprintf('[%s EVENT] SKIPPED (event %d) due to full occupancy\n',datetime('now'),i_event);
                        else
                            fprintf('[%s EVENT] SKIPPED (event %d) due to violating operationg hours\n',datetime('now'),i_event);
                        end
                    end
                end
            end
        end
    end
    
    % update agg
    keys = single.keys();
    if ~isempty(keys)
        for ev = keys
            if contains(ev{1},'EV')
                if  k <= single(ev{1}).time.end % is charging duration
                    TOU = interp1(0:0.25:24-0.25,par.TOU,k,'nearest');
                    if length(single(ev{1}).powers) > 1
                        power = interp1(linspace(single(ev{1}).time.start, ...
                                            single(ev{1}).time.end,...
                                            length(single(ev{1}).powers)), ...
                                            single(ev{1}).powers, k);
                    elseif length(single(ev{1}).powers) == 1
                        power = single(ev{1}).powers;
                    end
                    if power == opt.prb.station.pow_max % hyperthetically when power is max power it's the  uncontrol charging
                        sim_single.profit_charging_uc(i_k) = sim_single.profit_charging_uc(i_k) + par.Ts * power * (single(ev{1}).price - TOU);
                    else % flexible charging
                        sim_single.profit_charging_c(i_k) = sim_single.profit_charging_c(i_k) + par.Ts * power * (single(ev{1}).price - TOU);
                    end
                    sim_single.power(i_k) = sim_single.power(i_k) + power;
%                     sim.profit_charging(i_k) = sim.profit_charging(i_k) + par.Ts * power * (station(ev{1}).price - TOU);
                    sim_single.occ.charging(i_k) = sim_single.occ.charging(i_k) + 1;
                else % is overstaying
                    if k <= single(ev{1}).time.leave 
                        sim_single.profit_overstay(i_k) = sim_single.profit_overstay(i_k) + par.Ts * single(ev{1}).tariff.overstay;
                        sim_single.occ.overstay(i_k) = sim_single.occ.overstay(i_k) + 1;
%                         sim.overstay_duration(i_k) = sim.overstay_duration(i_k) + par.Ts;
                    else
                        single.remove(ev{1});
                        single('num_occupied_pole') = single('num_occupied_pole') - 1;
                        single('num_empty_pole') = single('num_empty_pole') + 1;
                    end 
                end
            elseif contains(ev{1},'occ')
                sim_single.occ.total(i_k) = single('num_occupied_pole');
                sim_single.occ.empty(i_k) = single('num_empty_pole');
            end
        end
    end
end

sim_single.par = par;

%% Simulation -- station wide
t = par.sim.starttime:par.Ts:par.sim.endtime; i_k = 0; i_event = 0;
sim_station = init_sim(t); % simulation result
station = containers.Map; % station monitor
station('num_occupied_pole') = 0; 
station('FLEX_list') = [];
station('ASAP_list') = [];
station('num_empty_pole') = par.station.num_poles;
station('D_init') = 0;
station('pow_cap') = par.station.num_poles * 6.6 ; % this value is arbitrary for now
station('cost_dc') = 200; % this value is arbitrary for now
sim_station.events = events;

for k = par.sim.starttime:par.Ts:par.sim.endtime
    i_k = i_k + 1;
    no_event_counter = 1;
    % check visit
    if i_event <= length(events.time)
        if any(round(events.time/par.Ts)*par.Ts == k)
            inds_events_k = find(round(events.time/par.Ts)*par.Ts == k);
            no_event_counter = 1; % reset counter
            
            for j = 1:length(inds_events_k)
                i_event = i_event + 1; % number of investigated events
                if events.inp{inds_events_k(j)}.duration <= par.sim.endtime - k ...
                        && station('num_empty_pole') > 0
                   sim_station.tot_decision = sim_station.tot_decision + 1; % number of decisions
                   sim_station.events.triggered(i_event) = true; % this event is triggered
                   
                   set_glob_prb(init_prb(events.inp{inds_events_k(j)}));

                   % find optimal tariff
                   if isempty(station('FLEX_list')) && isempty(station('ASAP_list'))
                       opt = run_opt();
                   else
                       [station, opt] = run_opt_station(station, k);
                   end
                   sim_station.opts{i_event} = opt;
%                    station('D_init') = opt.peak_pow; % update demand charge
                   
                   % driver makes choice
                   if par.sim.isFixedSeed
                       rng(1);
                   end
                   rc = 0; % XXX always choose flex
                   if rc <= opt.prob.flex
                       opt.choice = 0; % charging with flexibility
                       opt.time.end = opt.time.end_flex;
                       opt.powers = opt.flex.powers; % TODO: need to adjust the power profile for existing users with right time index
                       opt.price = opt.tariff.flex;
                       if isempty(station('FLEX_list')) % record FLEX user
                           station('FLEX_list') = {['EV' num2str(sim_station.tot_decision)]};
                       else
                           station('FLEX_list') = [station('FLEX_list') {['EV' num2str(sim_station.tot_decision)]}];
                       end
                   elseif rc <= opt.prob.flex + opt.prob.asap
                       opt.choice = 1; % charging as soon as possible
                       opt.time.end = opt.time.end_asap;
                       opt.powers = opt.asap.powers;
                       opt.price = opt.tariff.asap;
                       if isempty(station('ASAP_list')) % record FLEX user
                           station('ASAP_list') = {['EV' num2str(sim_station.tot_decision)]};
                       else
                           station('ASAP_list') = [station('ASAP_list') {['EV' num2str(sim_station.tot_decision)]}];
                       end
                   else
                       opt.choice = 2; % leaving without charging
                   end
                   sim_station.choice_probs(i_event,:) = opt.v;
                   sim_station.choice(i_event) = opt.choice;
                   sim_station.control(i_event,:) = opt.z(1:3);
                   if par.VIS_DETAIL
                    fprintf('[%s EVENT] time = %.2f, CHOICE = %s\n',datetime('now'),k,par.dcm.choices{opt.choice+1});
                   end

                   % if the driver chooses to charge EV
                   if opt.choice <= 1
                       [opt.time.leave, duration] = get_rand_os_duration(opt);
                       try
                           sim_station.overstay_duration(i_k) = sim_station.overstay_duration(i_k) + duration;
                           if duration == 32
                               disp('--- duration is boomed ----')
                           end
                       catch
                           a = 1;
                       end
                       sim_station.num_service(i_k) = sim_station.num_service(i_k) + 1;
                       station('num_occupied_pole') = station('num_occupied_pole') + 1;
                       station('num_empty_pole') = station('num_empty_pole') - 1;
                       opt.power_traj_actual = [];
                       station(['EV' num2str(sim_station.tot_decision)]) = opt;
                   end 
                else
                    if par.VIS_DETAIL
                        if station('num_empty_pole') == 0
                            fprintf('[%s EVENT] SKIPPED (event %d) due to full occupancy\n',datetime('now'),i_event);
                        else
                            fprintf('[%s EVENT] SKIPPED (event %d) due to violating operationg hours\n',datetime('now'),i_event);
                        end
                    end
                end
            end
        else
            no_event_counter = no_event_counter + 1;
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
                    power = station(ev{1}).powers(no_event_counter);
%                     power = station(ev{1}).powers(1);
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
                        sim_station.profit_charging_uc(i_k) = sim_station.profit_charging_uc(i_k) + par.Ts * power * (station(ev{1}).price - TOU);
                    else % flexible charging
                        sim_station.profit_charging_c(i_k) = sim_station.profit_charging_c(i_k) + par.Ts * power * (station(ev{1}).price - TOU);
                    end
                    
                    sim_station.power(i_k) = sim_station.power(i_k) + power;
%                     if station('D_init') < sim_station.power(i_k)
%                         % update demand charge
%                         station('D_init') = sim_station.power(i_k);
%                     end
                    
%                     sim.profit_charging(i_k) = sim.profit_charging(i_k) + par.Ts * power * (station(ev{1}).price - TOU);
                    sim_station.occ.charging(i_k) = sim_station.occ.charging(i_k) + 1;
                else % is overstaying
                    if k <= station(ev{1}).time.leave 
                        sim_station.profit_overstay(i_k) = sim_station.profit_overstay(i_k) + par.Ts * station(ev{1}).tariff.overstay;
                        sim_station.occ.overstay(i_k) = sim_station.occ.overstay(i_k) + 1;
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
                sim_station.occ.total(i_k) = station('num_occupied_pole');
                sim_station.occ.empty(i_k) = station('num_empty_pole');
            end
        end
    end
    if station('D_init') < sim_station.power(i_k)
        % update demand charge
        station('D_init') = sim_station.power(i_k);
    end
end

sim_station.par = par;


%% Baseline
t = par.sim.starttime:par.Ts:par.sim.endtime; i_k = 0; i_event = 0;
sim_base = init_sim(t); % simulation result
station = containers.Map; % station monitor
station('num_occupied_pole') = 0; 
station('num_empty_pole') = par.station.num_poles;

for k = par.sim.starttime:par.Ts:par.sim.endtime
    i_k = i_k + 1;
    % check visit
    if i_event <= length(events.time)
        if any(round(events.time/par.Ts)*par.Ts == k)
            inds_events_k = find(round(events.time/par.Ts)*par.Ts == k);
            for j = 1:length(inds_events_k)
                i_event = i_event + 1; % number of investigated events

                event = events.inp{inds_events_k(j)};
                if event.duration <= par.sim.endtime - k ...
                        && station('num_empty_pole') > 0
                   sim_base.tot_decision = sim_base.tot_decision + 1; % number of decisions
                   sim_base.events.triggered(i_event) = true; % this event is triggered

                   prb = set_glob_prb(init_prb(event));

                   opt.prb = prb;
                   opt.power = (event.SOC_need-event.SOC_init)*event.batt_cap/event.duration;
                   opt.time.start = k;
                   opt.time.end = k + event.duration;
                   opt.time.leave = k + event.duration + event.overstay_duration;
                   opt.tariff.overstay = par.base.tariff.overstay;

                   sim_base.num_service(i_k) = sim_base.num_service(i_k) + 1;
                   sim_base.overstay_duration(i_k) = sim_base.overstay_duration(i_k) + event.overstay_duration;
                   station('num_occupied_pole') = station('num_occupied_pole') + 1;
                   station('num_empty_pole') = station('num_empty_pole') - 1;
                   station(['EV' num2str(sim_base.tot_decision)]) = opt;
                else
                    if par.VIS_DETAIL
                        if station('num_empty_pole') == 0
                            fprintf('[%s EVENT] SKIPPED (event %d) due to full occupancy\n',datetime('now'),i_event);
                        else
                            fprintf('[%s EVENT] SKIPPED (event %d) due to violating operationg hours\n',datetime('now'),i_event);
                        end
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
                if  k <= station(ev{1}).time.end % is charging duration
                    TOU = interp1(0:0.25:24-0.25,par.TOU,k,'nearest');
                    power = station(ev{1}).power;
                    sim_base.power(i_k) = sim_base.power(i_k) + power;
                    sim_base.occ.charging(i_k) = sim_base.occ.charging(i_k) + 1;
                    sim_base.profit_charging_uc(i_k) = sim_base.profit_charging_uc(i_k) + par.Ts * power * (0.55*TOU); % the charging fee is twice as high as TOU
%                     if k >= station(ev{1}).time.start + 4
%                         sim.profit(i_k) = sim.profit(i_k) + par.Ts * station(ev{1}).tariff.overstay;
%                     end
                else % is overstaying
                    if k <= station(ev{1}).time.leave 
                        sim_base.occ.overstay(i_k) = sim_base.occ.overstay(i_k) + 1;
                        sim_base.profit_overstay(i_k) = sim_base.profit_overstay(i_k) + par.Ts * station(ev{1}).tariff.overstay;
%                         sim.overstay_duration(i_k) = sim.overstay_duration(i_k) + par.Ts;
                    else
                        station.remove(ev{1});
                        station('num_occupied_pole') = station('num_occupied_pole') - 1;
                        station('num_empty_pole') = station('num_empty_pole') + 1;
                    end 
                end
            elseif contains(ev{1},'occ')
                sim_base.occ.total(i_k) = station('num_occupied_pole');
                sim_base.occ.empty(i_k) = station('num_empty_pole');
            end
        end
    end
end


%% Visualize
viz = false;
if viz
    i=1;sim_station = sim_results_v2{i};sim_single=sim_results{i};sim_base=sim_results_base{i};
    plot(sim_station.t, sim_station.power, 'linewidth',1.5); hold on;
    plot(sim_single.t, sim_single.power, 'linewidth',1.5);   
    plot(sim_base.t, sim_base.power, 'linewidth',1.5); 
    plot(sim_station.t, 50*ones(1,length(sim_station.t)), '--', 'linewidth',1.5);hold off; 
    set(gca, 'fontsize',15); grid on; 
    legend('Station','Single','Baseline','capacity'); 
    xlabel('Time of the day [hr]'); ylabel('Power [kW]');
end

%% return state 1: max(single, station) = single, 0: max(single, station) = station
if max(sim_single.power) < max(sim_station.power)
    stop = 1;
end
varargout{1} = max(sim_single.power) > max(sim_station.power);
end