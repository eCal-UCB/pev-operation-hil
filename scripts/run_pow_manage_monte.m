% function varargout = run_pow_manage_monte(varargin)

% if nargin == 0
    par = set_glob_par(init_params());
% elseif nargin == 1
%     par = set_glob_par(varargin{1});
% end

% if par.sim.isFixedEventSequence
%     s = 'fixed';
% else
%     s = 'rand';
% end

num_sim = par.monte.num_sims; % simulation numbers

if par.sim.isFixedSeed
    seed_list = linspace(1, num_sim, num_sim);
%     seed_list = linspace(21,50,30);
end

failure_counter = 0;
n = 1;

while n <= num_sim
     
    try
        if par.sim.isFixedSeed
            seed_val = seed_list(n);
            events = gen_events_one_day(par, seed_val);
        else
            events = gen_events_one_day(par);
        end

        disp('Run power management for three settings.')
        management_state = run_power_management(par, events);
        failure_counter = failure_counter + management_state;

    catch e
        disp(e.message)
        warning('Error')
    end
    n = n + 1;
end

disp(['Out of ', num2str(num_sim), ' simulations, ', ...
    num2str(failure_counter), ' succeed in addressing peak power.'])