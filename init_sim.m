function sim = init_sim(t)
par = get_glob_par();

sim.power = zeros(size(t));
sim.profit = zeros(size(t));
sim.profit_charging_uc = zeros(size(t)); % record revenue from providing charging service uncontrolled charging
sim.profit_charging_c = zeros(size(t)); % record revenue from providing charging service flexible charging
sim.profit_overstay = zeros(size(t)); % record penalty from overstay
sim.occ.total  = zeros(size(t));
sim.occ.empty  = zeros(size(t));
sim.occ.charging = zeros(size(t));
sim.occ.overstay = zeros(size(t));
sim.overstay_duration = zeros(size(t));
sim.choice = nan*ones(par.sim.num_events,1);
sim.choice_probs = nan*ones(par.sim.num_events,3);
sim.num_service = zeros(size(t));
sim.tot_decision = 0;
sim.tot_num_service = 0;
sim.t = t;
sim.control = zeros(length(t),3);
sim.opts = cell(size(t));
end