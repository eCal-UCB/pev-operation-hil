function par = init_params()
% initialize parameters
par.operating_hours = [];                   % operating hours
par.station.num_poles = [];                 % number of charging poles
par.TOU = [];                               % TOU price
par.dcm.charging_flex.params = [];  % DCM parameters for choice 1 -- charging with flexibility 
par.dcm.charging_asap.params = [];  % DCM parameters for choice 2 -- charging as soon as possible
par.dcm.leaving.params = [];             % DCM parameters for choice 3 -- leaving without charging
par.pdf.visiting_rate = [];
end
