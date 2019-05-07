function events = gen_events_one_day()
% randomly generates a sequence of events
% TODO: generate events in order of ascending time.

act_data = readtable('real_act_data.csv');

num_events = 30;
% rng(1)
event_idx = sort(randi([1 height(act_data)], 1, num_events));

events.inp = cell(num_events,1);
events.time = zeros(num_events,1);

% test data -- to be removed
% test_times = [0:0.5:num_events] + 12;


for i = 1:num_events
    n = event_idx(i); % specify event index 
    if act_data{n, 6} < 0.3
        continue
    end
    event.time = act_data{n, 2}; 
    event.SOC_init = act_data{n, 3};
    event.SOC_need = act_data{n, 4}; % add infeasible scenario
    event.batt_cap = act_data{n, 5};
    event.duration = act_data{n, 6}; % hours
    event.overstay_duration = act_data{n, 7};
    event.pow_max = act_data{n, 8};
    event.pow_min = 0;
    
    events.inp{i} = event;
    events.time(i) = event.time;
end

end