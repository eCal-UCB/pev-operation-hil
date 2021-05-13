%% calculate actual power trajectory
power_profile = zeros(size(sim_station.t));
for i = 1:length(sim_station.opts)
    try
        opt = sim_station.opts{1,i};
        start_time_idx = opt.time.start/0.25;
        end_time_idx = opt.time.end/0.25;
        power_profile(start_time_idx:end_time_idx-1) = power_profile(start_time_idx:end_time_idx-1) + opt.power_traj_actual;
    catch
        continue
    end
end
%% plot energy
figure(10)
plot(cumsum(optSol.E_DA(1:24)),'LineWidth',2)
hold on

try
    plot(0:0.25:24-0.25, cumsum(power_profile(1:96)))
catch
    plot(0:0.25:24-0.25, cumsum(sim_station.power(1:96)))
end
hold off
ylabel('kWh','FontSize', 16)
xlabel('Time of Day','FontSize', 16)
legend('E DA', 'Sim Power', 'FontSize', 16)

figure(11)
plot(cumsum(optSol.E_DA(1:24)),'LineWidth',2)
hold on

try
    plot(cumsum(mean(reshape(power_profile(1:96), 4, []))),'LineWidth',2)
catch
    plot(cumsum(mean(reshape(sim_station.power(1:96), 4, []))))
end
hold off
ylabel('kWh','FontSize', 16)
xlabel('Time of Day','FontSize', 16)
legend('E DA', 'Sim Power', 'FontSize', 16)

%%
power = 0;
for i = 1:length(sim_station.opts)
    try
        n = (sim_station.opts{1,i}.var_dim_constant-1)/2;
        power = power + sum(sim_station.opts{1,i}.x(n+2:end));
    catch
        continue
    end
end
power = power * 0.25;


