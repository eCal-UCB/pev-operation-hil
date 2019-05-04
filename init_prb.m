function prb = init_prb(varargin)
par = get_glob_par();

if nargin == 0
    % user input
    prb.user.time     = 7.25;
    prb.user.SOC_init = 0.2;
    prb.user.SOC_need = 0.5;
    prb.user.batt_cap = 100; % kwh
    prb.user.duration = 8; % hrs
elseif nargin == 1
    inp = varargin{1};
    prb.user.time     = inp.time;
    prb.user.SOC_init = inp.SOC_init;
    prb.user.SOC_need = inp.SOC_need;
    prb.user.batt_cap = inp.batt_cap;
    prb.user.duration = inp.duration;
end

% problem specifications
prb.N_flex = prb.user.duration/par.Ts;             % charging duration that is not charged, hour
prb.N_asap = ceil((prb.user.SOC_need-prb.user.SOC_init)...
             *prb.user.batt_cap/par.station.pow_max/par.eff/par.Ts);
prb.TOU = interp1(0:23,par.TOU,...
        prb.user.time:par.Ts:prb.user.time+prb.user.duration-par.Ts,... 
        'nearest')';                               % TOU price
end