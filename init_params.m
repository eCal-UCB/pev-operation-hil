function par = init_params(varargin)

% sensitivity analysis
% par.sens_analysis.num_poles = 2:2:16;
par.sens_analysis.num_poles = [2,3];

% monte carlo
par.monte.num_sims = 10;

% each one day simulation
par.sim.starttime = 7;
par.sim.endtime = 22;
par.sim.isFixedEventSequence = false;
par.sim.num_events = 40;
par.Ts = 0.25; % timestep, hour -- must decompose 1

% baseline parameters
par.base.tariff.overstay = 1.0;

% TOU
par.TOU = [0.217*ones(1,34) ...    % 0-8.5
           0.244*ones(1,48-34) ...    % 8.5-12
           0.268*ones(1,72-48) ...      % 12-16
           0.244*ones(1,86-72) ...    % 16-21.5
           0.217*ones(1,96-86)];      % 22-24

% charging station config
par.station.num_poles = 7;                 % number of charging poles
par.eff = 0.89;                             % power efficiency

% dcm params
par.dcm.choices = [{'charging with flexibility'},{'charging asap'},{'leaving without charging'}];

% pdfs
par.pdf.visit = [0.1*ones(1,7) ...    % 0-7
                 0.3*ones(1,5) ...    % 7-12
                 0.2*ones(1,2) ...    % 12-14
                 0.2*ones(1,2) ...    % 14-16
                 0.2*ones(1,6) ...    % 16-22
                 0.001*ones(1,2)];    % 22-24

% regularization params
par.lambda.x = 10;
par.lambda.z_c = 10;
par.lambda.z_uc = 0.1;
par.lambda.h_c = 0.01; % TODO: should be average overstay penalty in real data, should move to par
par.lambda.h_uc = 0.1; % TODO: should be average overstay penalty in real data, should move to par
par.mu = 1e4;
par.soft_v_eta = 1e-3; % softening equality constraint for v; to avoid numerical error
par.opt.eps = 1e-4;
end
