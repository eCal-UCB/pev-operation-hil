function events = gen_events_one_day()
% randomly generates a sequence of events
% TODO: generate events in order of ascending time.

num_events = 10;
events.inp = cell(num_events,1);
events.time = zeros(num_events,1);

% test data -- to be removed
test_times = [1:num_events] + 7;


for n = 1:num_events
    % TODO: generate event structure based on the real data set
    event.time = test_times(n); 
    event.SOC_init = 0.2;
    event.SOC_need = 0.4; % add infeasible scenario
    event.batt_cap = 24;
    event.duration = 3; % hours
    event.overstay_duration = 2.3;
    event.pow_max = 3.2;
    event.pow_min = 0;
    
    events.inp{n} = event;
    events.time(n) = event.time;
end

end