function par = init_params()
% initialize parameters

par.station.num_poles = 10;                 % number of charging poles
par.station.isAllOccupied = true;
par.station.pow_min = 0;
par.station.pow_max = 7.2;                  % kw

par.TOU = [];                               % TOU price

par.dcm.charging_flex.params = [];          % DCM parameters for choice 1 -- charging with flexibility 
par.dcm.charging_asap.params = [];          % DCM parameters for choice 2 -- charging as soon as possible
par.dcm.leaving.params = [];                % DCM parameters for choice 3 -- leaving without charging

par.pdf.visiting_rate = [];

par.user.SOC_init = [];
par.user.batt_cap = [];
par.user.SOC_need = [];
par.user.duration = [];

par.Ts = 1; % timestep, hour
par.N_flex = par.user.duration;             % charging duration that is not charged, hour
par.N_asap = ceil((par.user.SOC_need-par.user.SOC_init)...
             *par.user.batt_cap/par.station.max_pow);
par.lambda.h_c = 10;
par.lambda.h_uc = 10;
end
