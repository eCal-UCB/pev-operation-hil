function prb = init_prb(varargin)
par = get_glob_par();

if nargin == 0
    % user input
    prb.user.time     = 14.25;
    prb.user.SOC_init = 0.3;
    prb.user.SOC_need = 0.5;
    prb.user.batt_cap = 80; % kwh
    prb.user.duration = 8; % hrs
    prb.user.overstay_duration = 1;
    prb.station.pow_max = 7.2;
    prb.station.pow_min = 0;
elseif nargin == 1
    event = varargin{1};
    prb.user.time     = round(event.time/par.Ts)*par.Ts;
    prb.user.SOC_init = event.SOC_init;
    prb.user.SOC_need = event.SOC_need;
    prb.user.batt_cap = event.batt_cap;
    prb.user.duration = round(event.duration/par.Ts)*par.Ts;
    prb.user.overstay_duration = round(event.overstay_duration/par.Ts)*par.Ts;
    prb.station.pow_max = event.pow_max;
    prb.station.pow_min = event.pow_min;
end

% problem specifications
prb.N_flex = prb.user.duration/par.Ts;             % charging duration that is not charged, hour
prb.N_asap = ceil((prb.user.SOC_need-prb.user.SOC_init)...
             *prb.user.batt_cap/prb.station.pow_max/par.eff/par.Ts);
prb.TOU = interp1(0:0.25:24-0.25,par.TOU,...
        prb.user.time:par.Ts:prb.user.time+prb.user.duration-par.Ts,... 
        'nearest')';                               % TOU price
end