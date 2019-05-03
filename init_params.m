function par = init_params()
% initialize parameters

% user input
par.user.SOC_init = 0.2;
par.user.SOC_need = 0.5;
par.user.batt_cap = 100; % kwh
par.user.duration = 8; % hrs

% charging station config
par.station.num_poles = 10;                 % number of charging poles
par.station.isAllOccupied = true;
par.station.pow_min = 0;
par.station.pow_max = 7.2;                  % kw
par.eff = 0.89;                             % power efficiency

% dcm params
par.dcm.charging_flex.params = [-0.01 0 0 1]';          % DCM parameters for choice 1 -- charging with flexibility 
par.dcm.charging_asap.params = [0 -0.01 0 1.5]';          % DCM parameters for choice 2 -- charging as soon as possible
par.dcm.leaving.params       = [0 0 0.01 0]';           % DCM parameters for choice 3 -- leaving without charging
par.THETA = [par.dcm.charging_flex.params';
             par.dcm.charging_asap.params';
             par.dcm.leaving.params'];

% pdfs
par.pdf.visiting_rate = [];

% problem specifications
par.Ts = 0.2; % timestep, hour -- must decompose 1
par.N_flex = par.user.duration/par.Ts;             % charging duration that is not charged, hour
par.N_asap = ceil((par.user.SOC_need-par.user.SOC_init)...
             *par.user.batt_cap/par.station.pow_max/par.eff/par.Ts);
par.opt.eps = 1e-4;

% regularization params
par.lambda.h_c = 3;
par.lambda.h_uc = 3;
par.mu = 1e4;
par.soft_v_eta = 0.001; % softening equality constraint for v; to avoid numerical error

% TOU price
par.TOU = 1./cumsum(ones(par.user.duration/par.Ts,1));                               % TOU price
end
