function par = init_params()
% initialize parameters
par.operating_hours = 10;                   % operating hours
par.station.num_poles = 10;                 % number of charging poles
par.TOU = [];                               % TOU price
par.dcm.charging_flex.params = [-0.1 0.2];  % DCM parameters for choice 1 -- charging with flexibility 
par.dcm.charging_asap.params = [-0.1 0.2];  % DCM parameters for choice 2 -- charging as soon as possible
par.dcm.leaving.params = [0 0];             % DCM parameters for choice 3 -- leaving without charging
end
