function sim = run_sim_one_day_baseline(sim_c)
% input: a one-day simulation result with controller
fprintf('[%s SIM] start running simulation with baseline (without controller)\n',datetime('now'));
events = sim_c.events; par = sim_c.par;
t = par.sim.starttime:par.Ts:par.sim.endtime; i_k = 0; i_event = 0;
sim = init_sim(t); % simulation result
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
                   sim.tot_decision = sim.tot_decision + 1; % number of decisions
                   sim.events.triggered(i_event) = true; % this event is triggered

                   prb = set_glob_prb(init_prb(event));

                   opt.prb = prb;
                   opt.power = (event.SOC_need-event.SOC_init)*event.batt_cap/event.duration;
                   opt.time.start = k;
                   opt.time.end = k + event.duration;
                   opt.time.leave = k + event.duration + event.overstay_duration;
                   opt.tariff.overstay = par.base.tariff.overstay;

                   sim.num_service(i_k) = sim.num_service(i_k) + 1;
                   station('num_occupied_pole') = station('num_occupied_pole') + 1;
                   station('num_empty_pole') = station('num_empty_pole') - 1;
                   station(['EV' num2str(sim.tot_decision)]) = opt;
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
                if  k <= station(ev{1}).time.end % is charging duration
                    TOU = interp1(0:0.25:24-0.25,par.TOU,k,'nearest');
                    power = station(ev{1}).power;
                    sim.power(i_k) = sim.power(i_k) + power;
                    sim.occ.charging(i_k) = sim.occ.charging(i_k) + 1;
                    sim.profit_charging_uc(i_k) = sim.profit_charging_uc(i_k) + par.Ts * power * (TOU); % the charging fee is twice as high as TOU
%                     if k >= station(ev{1}).time.start + 4
%                         sim.profit(i_k) = sim.profit(i_k) + par.Ts * station(ev{1}).tariff.overstay;
%                     end
                else % is overstaying
                    if k <= station(ev{1}).time.leave 
                        sim.occ.overstay(i_k) = sim.occ.overstay(i_k) + 1;
                        sim.profit_overstay(i_k) = sim.profit_overstay(i_k) + par.Ts * station(ev{1}).tariff.overstay;
                        sim.overstay_duration(i_k) = sim.overstay_duration(i_k) + par.Ts;
                    else
                        station.remove(ev{1});
                        station('num_occupied_pole') = station('num_occupied_pole') - 1;
                        station('num_empty_pole') = station('num_empty_pole') + 1;
                    end 
                end
            elseif contains(ev{1},'occ')
                sim.occ.total(i_k) = station('num_occupied_pole');
                sim.occ.empty(i_k) = station('num_empty_pole');
            end
        end
    end
end
end