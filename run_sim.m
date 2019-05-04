% This script is to simulate EV charging station operations where the
% charging tariff is determined real-time with taking account into EV
% drivers' behaviors. The overall objective of the tariff control is to
% minimize the operation cost of the system operater. 
%
% There are three choices that each EV driver can make at arrival:
% (i) charging with flexibility, (ii) charging as soon as possible, and
% (iii) leaving without charging.
%
% At each arrival of EV driver, user specific parameters, e.g., battery
% capacity, desired parking durations, initial SOC, and needed SOC level
% (for next mobility demand) are randomly sampled from an empirical
% probability distribution function that is generated with TELD dataset 
% [X].
%
% THIS WORK IS A PART OF EE227C COURSE PROJECT AT UC BERKELEY.
% last modified, May 2019.

%% Initialization
disp('[ INIT] initializing...');
par = set_glob_par(init_params());
disp('[ INIT] DONE');

%% Simulation

for k = par.sim.starttime:par.Ts:par.sim.endtime
    
    % random visit
    r = rand;
    if r < interp1(0:23,par.pdf.visit,k)
       % TODO: sample values from pdf
       inp.user.time = k;
       inp.SOC_init = 0.2;
       inp.SOC_need = 0.5;
       inp.batt_cap = 24;
       inp.duration = 6;
       set_glob_prb(init_prb(inp));
       
       % run opt
       opt = run_opt();
    end
end
%% Visualization