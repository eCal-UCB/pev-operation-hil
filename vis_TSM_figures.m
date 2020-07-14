% This script visualizes figures that are included in 
% "inducing human bheavior to maximize operatino performance at PEV charging stations"
% submitted to IEEE Transaction on Smart Grid.
%
% The figures that are visualized through this scripts include:
% 1. power management (station-wide vs single)
% 2. temporal electricity price changes (station-wide vs baseline)
% 3. sensitivity analysis to decompose profit

%% 1. Power management (station vs single)
% choose the day where:
% - sigle charger violates power capacity constraint, while station-wide
% satisfies
% - gap in max power is huge enough
load('monte-sim-results/07_12_20_10_44_monte_eps100_rand_seq_poles15.mat');
par = sim_results{1}.par;
N = length(sim_results);
power_cap = par.station.num_poles * 6.6;
violated_inds = [];
for n = 1 : N
    if max(sim_results{n}.power) >= power_cap ...
    && max(sim_results_v2{n}.power) <= power_cap
        violated_inds = [violated_inds, n];
    end
end

max_gap = 0;
max_ind = 0;
for n = violated_inds
    if max_gap <= max(sim_results_v2{n}.power) - max(sim_results{n}.power)
        max_gap = max(sim_results_v2{n}.power) - max(sim_results{n}.power);
        max_ind = violated_inds(n);
    end
end

figure;
plot(sim_results_v2{n}.t, sim_results_v2{n}.power,'b','linewidth',1.5); hold on;
plot(sim_results{n}.t, sim_results{n}.power,'r','linewidth',1.5); 
plot(sim_results_v2{n}.t, power_cap*ones(1,length(sim_results_v2{n}.power)),'k--', 'linewidth',1.5); hold off;
set(gca, 'fontsize', 15); grid on;
title('Power Management');
xlabel('hour of the day');
ylabel('power (kW)');
legend('station','single','capacity');


%% 2. Temporal electricity price changes (station-wide vs baseline)
load('monte-sim-results/06_26_20_18_21_monte_eps200_rand_seq_poles8.mat'); % this will overwrite workspace
day = 1;

N = length(sim_results_v2{day}.opts);
event_times = [];
flex_prices = [];
asap_prices = [];
overstay_prices = [];
TOU = []; % TODO: interpolate
for i = 1 : N
    opt = sim_results_v2{day}.opts{i};
    if ~isempty(opt)
        event_times = [event_times, opt.time.start];
        flex_prices = [flex_prices, opt.z(1)];
        asap_prices = [asap_prices, opt.z(2)];
        overstay_prices = [overstay_prices, opt.z(3)];
        TOU = [TOU, interp1(0:0.25:24-0.25, par.TOU, opt.time.start)];
    end
end

figure;
colororder({'b','m'})

yyaxis left
p1 = plot(event_times, flex_prices, 'linewidth', 1.5); hold on;
p2 = plot(event_times, asap_prices, 'linewidth', 1.5);
p4 = plot(event_times, overstay_prices, 'linewidth', 1.5); hold off;
ylabel('price ($)');

yyaxis right
p3 = plot(event_times, TOU, 'linewidth', 1.5); hold off;
set(gca, 'fontsize', 15); grid on;
title('Charging Price');
xlabel('hour of the day');
ylabel('price ($)');
legend('flex','asap','overstay','TOU');


%% 3. Stated duration vs energy requested vs overstay penalty
load('monte-sim-results/06_26_20_18_21_monte_eps200_rand_seq_poles8.mat'); % this will overwrite workspace

durations = [];
energy_requested = [];
overstay_penalty = [];
for day = 1 : length(sim_results)
    N = length(sim_results{day}.events.inp);
    for i = 1 : N
        if ~isempty(sim_results{day}.opts{i})
            inp = sim_results{day}.events.inp{i};
            durations = [durations, inp.duration];
            energy_requested = [energy_requested, (inp.SOC_need - inp.SOC_init) * inp.batt_cap];
            overstay_penalty = [overstay_penalty, sim_results{day}.opts{i}.z(3)];
        end
    end
end

figure(1); scatter(energy_requested,durations, 50, overstay_penalty, 'filled');
grid on; 
colormap cool

cb = colorbar();
set(gca,'fontsize',15);
xlabel('energy requested (kWh)');
ylabel('stated duration (hr)');
legend('overstay penalty ($/hr)');