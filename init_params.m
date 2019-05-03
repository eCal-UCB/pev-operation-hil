function par = init_params()
% initialize parameters

% charging station config
par.station.num_poles = 10;                 % number of charging poles
par.station.isAllOccupied = true;
par.station.pow_min = 0;
par.station.pow_max = 7.2;                  % kw

% dcm params
par.dcm.charging_flex.params = [];          % DCM parameters for choice 1 -- charging with flexibility 
par.dcm.charging_asap.params = [];          % DCM parameters for choice 2 -- charging as soon as possible
par.dcm.leaving.params = [];                % DCM parameters for choice 3 -- leaving without charging

% pdfs
par.pdf.visiting_rate = [];

% user input
par.user.SOC_init = 0.2;
par.user.SOC_need = 0.5;
par.user.batt_cap = 24; % kwh
par.user.duration = 8; % hrs

% problem specifications
par.Ts = 1; % timestep, hour
par.N_flex = par.user.duration;             % charging duration that is not charged, hour
par.N_asap = ceil((par.user.SOC_need-par.user.SOC_init)...
             *par.user.batt_cap/par.station.max_pow);

% regularization params
par.lambda.h_c = 10;
par.lambda.h_uc = 10;

% TOU price
par.TOU = ones(par.N_flex,1);                               % TOU price
end
