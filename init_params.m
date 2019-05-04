function par = init_params(varargin)

% simulation parameters
par.sim.starttime = 7;
par.sim.endtime = 22;
par.Ts = 0.25; % timestep, hour -- must decompose 1

% TOU
par.TOU = [0.2*ones(1,7) ...    % 0-7
           0.5*ones(1,5) ...    % 7-12
           1*ones(1,2) ...      % 12-14
           2*ones(1,2) ...      % 14-16
           0.7*ones(1,6) ...    % 16-22
           0.2*ones(1,2)];      % 22-24

% charging station config
par.station.num_poles = 10;                 % number of charging poles
par.station.isAllOccupied = true;
par.station.pow_min = 0;
par.station.pow_max = 7.2;                  % kw
par.eff = 0.89;                             % power efficiency

% dcm params
par.dcm.choices = [{'charging with flexibility'},{'charging asap'},{'leaving without charging'}];
par.dcm.charging_flex.params = [-0.01 0 0 1]';          % DCM parameters for choice 1 -- charging with flexibility 
par.dcm.charging_asap.params = [0 -0.01 0 1.5]';          % DCM parameters for choice 2 -- charging as soon as possible
par.dcm.leaving.params       = [0 0 0.01 0]';           % DCM parameters for choice 3 -- leaving without charging
par.THETA = [par.dcm.charging_flex.params';
             par.dcm.charging_asap.params';
             par.dcm.leaving.params'];

% pdfs
par.pdf.visit = [0.1*ones(1,7) ...    % 0-7
                 0.5*ones(1,5) ...    % 7-12
                 0.2*ones(1,2) ...    % 12-14
                 0.2*ones(1,2) ...    % 14-16
                 0.2*ones(1,6) ...    % 16-22
                 0.001*ones(1,2)];    % 22-24

% regularization params
par.lambda.h_c = 3;
par.lambda.h_uc = 3;
par.mu = 1e4;
par.soft_v_eta = 0.001; % softening equality constraint for v; to avoid numerical error
par.opt.eps = 1e-4;
end
