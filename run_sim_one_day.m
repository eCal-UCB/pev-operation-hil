function varargout = run_sim_one_day(varargin)
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
% THIS WORK IS A PART OF EE227C COURSE PROJECT AT UC BERKELEY.
% last modified, May 2019.
%
% Contributors: Sangjae Bae, Teng Zeng, Bertrand Travacca.

% clear
%% Initialization
fprintf('[%s INIT] initializing...\n',datetime('now'));
if nargin == 0
    par = set_glob_par(init_params());
    events = gen_events_one_day();
elseif nargin == 1
    par = varargin{1};
    events = gen_events_one_day();
elseif nargin == 2
    par = varargin{1};
    events = varargin{2};
else
    error(sprintf('[%s ERROR] too many input arguments',datetime('now')));
end

fprintf('[%s INIT] DONE\n',datetime('now'));
%% Simulation
t = par.sim.starttime:par.Ts:par.sim.endtime;
i = 0; 
sim = init_sim(t); % simulation result
sim.events = events;
station = containers.Map; % station monitor
station('num_occupied_pole') = 0; 
ind_event = 1;

for k = par.sim.starttime:par.Ts:par.sim.endtime
    i = i+1;
    
    % check visit
    if ind_event <= length(events.time)
        if any(round(events.time/par.Ts)*par.Ts == k)
            inds_event = find(round(events.time/par.Ts)*par.Ts == k);
            for j = 1:length(inds_event)
                if events.inp{inds_event(j)}.duration < par.sim.endtime - k
                   set_glob_prb(init_prb(events.inp{inds_event(j)}));

                   % find optimal tariff
                   opt = run_opt();
                   sim.opts{i} = opt;

                   % driver makes choice
                   rc = rand;
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
                   sim.choice_probs(i,:) = opt.v;
                   sim.choice(i) = opt.choice;
                   sim.control(i,:) = opt.z(1:3);
                   fprintf('[%s EVENT] time = %.2f, CHOICE = %s\n',datetime('now'),k,par.dcm.choices{opt.choice+1});

                   % if the driver chooses to charge EV
                   if opt.choice <= 1
                       [opt.time.leave, duration] = get_rand_os_duration(opt);
                       sim.overstay_duration(i) = sim.overstay_duration(i) + duration;
                       station('num_occupied_pole') = station('num_occupied_pole') + 1;
                       station(['EV' num2str(sim.tot_visit)]) = opt;
                   end
                   sim.tot_visit = sim.tot_visit + 1;
                end
                ind_event = ind_event + 1;
            end
        end
    end
    
    % update agg
    keys = station.keys();
    if ~isempty(keys)
        for ev = keys
            if contains(ev{1},'EV')
                if  k <= station(ev{1}).time.end % is charging duration
                    TOU = interp1(0:23,par.TOU,k,'nearest');
                    if length(station(ev{1}).powers) > 1
                        power = interp1(linspace(station(ev{1}).time.start, ...
                                            station(ev{1}).time.end,...
                                            length(station(ev{1}).powers)), ...
                                            station(ev{1}).powers, k);
                    elseif length(station(ev{1}).powers) == 1
                        power = station(ev{1}).powers;
                    end
                    
                    sim.power(i) = sim.power(i) + power;
                    sim.profit(i) = sim.profit(i) + par.Ts * power * (station(ev{1}).price - TOU);
                    sim.charging(i) = sim.charging(i) + 1;
                else % is overstaying
                    if k <= station(ev{1}).time.leave 
                        sim.profit(i) = sim.profit(i) + par.Ts * station(ev{1}).tariff.overstay;
                        sim.overstay(i) = sim.overstay(i) + 1;
                    else
                        station.remove(ev{1});
                        station('num_occupied_pole') = station('num_occupied_pole') - 1;
                    end 
                end
            elseif contains(ev{1},'occ')
                sim.occ(i) = station('num_occupied_pole');
            end
        end
    end
end

varargout = {};
varargout{1} = sim;

%% Visualization
if nargout == 0
    options = vis_sim_one_day(); % options: display, temporals, choices
    options.temporals = true;
    vis_sim_one_day(sim,options);
end

end