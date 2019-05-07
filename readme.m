% <simulations>
% run_sim_monte.m   : run monte carlo simulations
% run_sim_one_day.m : run simulation of one day operation
% run_sim_one_event.m : run simulation of one control event
% get_rand_os_duration.m  : randomly sample an overstay duration
% init_sim          : initialize sim structure that stores one day
%                     operation simulation result
%
%
% <optimizations>
% run_opt           : solve optimization with block coordinate descent 
%                     (BCD) algorithm for a current EV
% argmin_z          : update zk in BCD algorithm
% argmin_x          : update xk in BCD algorithm
% argmin_v          : update vk in BCD algorithm
% init_prb          : initialize prb structure that stores event specific
%                     parameteres 
% init_par          : initialize par structure that stores constant
%                     (hyper) parameters
%
%
% <visualizations>
% vis_sim_monte.m   : visualize a result of monte carlo simulations results.
%                     Run this file to open a file brower and select a .mat 
%                     file in /monte-sim-results 
% vis_sim_one_day.m : visualize a simulation result of one day operation
% vis_sim_one_event.m : visualize a simulation result of one control event



