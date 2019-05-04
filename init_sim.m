function sim = init_sim(t)
sim.power = zeros(size(t));
sim.profit  = zeros(size(t));
sim.occ  = zeros(size(t));
sim.charging = zeros(size(t));
sim.overstay = zeros(size(t));
sim.overstay_duration = zeros(size(t));
sim.choice = nan*ones(size(t));
sim.choice_probs = zeros(length(t),3);
sim.tot_visit = 0;
sim.t = t;
end