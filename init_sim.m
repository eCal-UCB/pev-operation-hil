function sim = init_sim(t)
par = get_glob_par();

sim.power = zeros(size(t));
sim.profit  = zeros(size(t));
sim.occ  = zeros(size(t));
sim.charging = zeros(size(t));
sim.overstay = zeros(size(t));
sim.overstay_duration = zeros(size(t));
sim.choice = nan*ones(par.num_events,1);
sim.choice_probs = nan*ones(par.num_events,3);
sim.tot_visit = 0;
sim.t = t;
sim.control = zeros(length(t),3);
sim.opts = cell(size(t));
end