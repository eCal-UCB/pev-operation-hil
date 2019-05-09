function varargout = gen_events_one_day(varargin)
% -- scenario=0 for baseline information; num_events could be anything
% arbitrary
% -- scenario=1 for randomly generates a sequence of events; need to specify
% num_events

if nargin == 0
    par = get_glob_par();
elseif nargin == 1
    par = varargin{1};
else
    fprintf('[%s ERROR] invalid number of inputs',datetime('now'));
end

% initialize
if par.sim.isFixedEventSequence % with one fixed sequence of events
    act_data = readtable('real_act_data_1day.csv');
    num_events = height(act_data);
    event_idx = 1:num_events;
else % with multi random sequences of events
    num_events = par.sim.num_events; % by default. 
    act_data = readtable('real_act_data.csv');
    event_idx = sort(randi([1 height(act_data)], 1, num_events));
end

% build up events structure
events.inp = cell(num_events,1);
events.time = zeros(num_events,1);
events.triggered = false*ones(num_events,1); % triggered event flag
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

if nargout == 1
    varargout = {};
    varargout{1} = events;
else
    
end
end