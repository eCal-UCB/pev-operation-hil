function vis_sensitivity_analysis(varargin)
if nargin == 0
    [fname, fpath] = uigetfile;
    data = load(fullfile(fpath,fname));
    monte_results = data.monte_results; 
elseif nargin == 1
    monte_results = varargin{1};
end

%% Data process
num_monte = length(monte_results);
par = monte_results{1}.optimal_v2{1}.par;
num_sims = par.monte.num_sims;
len_overstay_duration = length(monte_results{1}.optimal_v2{1}.overstay_duration);

num_poles = par.sens_analysis.num_poles;
mean_profit_charging_uc_gap = zeros(num_monte,1);
mean_profit_charging_c_gap = zeros(num_monte,1);
mean_profit_overstay_gap = zeros(num_monte,1);
mean_profit_gap = zeros(num_monte,1);
% -- temporary for value comparison
profit_charging_uc_sum = zeros(num_monte,1);
profit_charging_c_sum = zeros(num_monte,1);
profit_charging_uc_sum_base = zeros(num_monte,1);
profit_overstay_sum = zeros(num_monte,1);
profit_overstay_sum_base = zeros(num_monte,1);
cost_demand_charge = zeros(num_monte,1);
cost_demand_charge_base = zeros(num_monte,1);
%
mean_overstay_gap = zeros(num_monte,1);
box_overstay_gap = zeros(num_sims, num_monte);
plot_service_tot = zeros(num_monte,1);
plot_service_tot_var = zeros(num_monte,1);
plot_service_tot_max = zeros(num_monte,1);
plot_service_tot_min = zeros(num_monte,1);
plot_service_tot_base = zeros(num_monte,1);
plot_service_tot_var_base = zeros(num_monte,1);
plot_service_tot_max_base = zeros(num_monte,1);
plot_service_tot_min_base = zeros(num_monte,1);
box_service_tot = zeros(num_sims, num_monte);
box_service_tot_base = zeros(num_sims, num_monte);
for i = 1:num_monte
    sim_results = monte_results{i}.optimal_v2; %%%%%%%%%%% TEMP: set to optimal_v2 for station wide optimization results
    sim_results_base = monte_results{i}.baseline;
    
    num_sims = sim_results{1}.par.monte.num_sims;
    
    % with controller
    overstay_tot = zeros(num_sims,1);
    overstay_mean = zeros(num_sims,1);
    overstay_penalty_mean = zeros(num_sims,1);
    profit_tot = zeros(num_sims,1);
    profit_charging_uc = zeros(num_sims,1);
    profit_charging_c = zeros(num_sims,1);
    profit_overstay = zeros(num_sims,1);
    service_tot = zeros(num_sims,1);
    pow_max = zeros(num_sims,1);
    for n = 1:num_sims
        if ~isempty(sim_results{n}.overstay_duration(sim_results{n}.overstay_duration~=0))
            overstay_mean(n) = mean(sim_results{n}.overstay_duration(sim_results{n}.overstay_duration~=0));
        else
            overstay_mean(n) = 0;
        end
        overstay_penalty_mean(n) = mean(sim_results{n}.control(sim_results{n}.control(:,3)~=0,3));
        overstay_tot(n) = sum(sim_results{n}.overstay_duration);
        profit_charging_uc(n) = sum(sim_results{n}.profit_charging_uc);
        profit_charging_c(n) = sum(sim_results{n}.profit_charging_c);
        profit_overstay(n) = sum(sim_results{n}.profit_overstay);
        profit_tot(n) = profit_charging_uc(n) + profit_charging_c(n) + profit_overstay(n);
        service_tot(n) = sum(sim_results{n}.num_service);
        pow_max(n) = max(sim_results{n}.power);
    end

    % without controller
    overstay_tot_base = zeros(num_sims,1);
    overstay_mean_base = zeros(num_sims,1);
    overstay_penalty_mean_base = zeros(num_sims,1);
    profit_tot_base = zeros(num_sims,1);
    profit_charging_uc_base = zeros(num_sims,1);
    profit_charging_c_base = zeros(num_sims,1);
    profit_overstay_base = zeros(num_sims,1);
    service_tot_base = zeros(num_sims,1);
    pow_max_base = zeros(num_sims,1);
    for n = 1:num_sims
        overstay_mean_base(n) = mean(sim_results_base{n}.overstay_duration(sim_results_base{n}.overstay_duration~=0));
        overstay_penalty_mean_base(n) = mean(sim_results_base{n}.control(sim_results_base{n}.control(:,3)~=0,3));
        overstay_tot_base(n) = sum(sim_results_base{n}.overstay_duration);
        profit_charging_uc_base(n) = sum(sim_results_base{n}.profit_charging_uc);
        profit_charging_c_base(n) = sum(sim_results_base{n}.profit_charging_c);
        profit_overstay_base(n) = sum(sim_results_base{n}.profit_overstay);
        profit_tot_base(n) = profit_charging_uc_base(n) + profit_charging_c_base(n) + profit_overstay_base(n);
        service_tot_base(n) = sum(sim_results_base{n}.num_service);
        pow_max_base(n) = max(sim_results_base{n}.power);
    end
    
    % compute the improvements
    box_overstay_gap(:,i) = (overstay_mean./overstay_mean_base-1)*100;
    mean_profit_gap(i) = (mean(profit_tot)/mean(profit_tot_base)-1)*100;
    mean_profit_charging_uc_gap(i) = (mean(profit_charging_uc)/mean(profit_charging_uc_base)-1)*100;
    mean_profit_charging_c_gap(i) = (mean(profit_charging_c)/mean(profit_charging_c_base)-1)*100;
    mean_profit_overstay_gap(i) = (mean(profit_overstay)/mean(profit_overstay_base)-1)*100;
    % -- temporary for value comparison
    profit_charging_uc_sum(i) = sum(profit_charging_uc);
    profit_charging_c_sum(i) = sum(profit_charging_c);
    profit_charging_uc_sum_base(i) = sum(profit_charging_uc_base);
    profit_overstay_sum(i) = sum(profit_overstay);
    profit_overstay_sum_base(i) = sum(profit_overstay_base);
    cost_demand_charge(i) = 18.86*max(pow_max);
    cost_demand_charge_base(i) = 18.86*max(pow_max_base);
    %
    mean_overstay_gap(i) = (mean(overstay_mean(~isnan(overstay_mean)))/mean(overstay_mean_base)-1)*100;
    plot_service_tot(i) = mean(service_tot)/par.sim.num_events;
    plot_service_tot_var(i) = var(service_tot/par.sim.num_events);
    plot_service_tot_max(i) = prctile(service_tot/par.sim.num_events,75);
    plot_service_tot_min(i) = prctile(service_tot/par.sim.num_events,25);
    plot_service_tot_base(i) = mean(service_tot_base/par.sim.num_events);
    plot_service_tot_var_base(i) = var(service_tot_base/par.sim.num_events);
    plot_service_tot_max_base(i) = prctile(service_tot_base/par.sim.num_events,75);
    plot_service_tot_min_base(i) = prctile(service_tot_base/par.sim.num_events,25);
    box_service_tot(:,i) = service_tot;
    box_service_tot_base(:,i) = service_tot_base;
end


%% Visualization
style = 2;
% style1: both bars on the same figure
if style == 1 % 09.21.2019 Teng: not change since not using. 09.20.2019, changed
    close;figure; 
    b1=bar(num_poles,mean_overstay_gap,'facecolor',[0.4660 0.6740 0.1880]); hold on;
    b2=bar(num_poles,mean_profit_charging_gap,'FaceColor',[0,0.7,0.9]); hold off;
%     b3=bar(num_poles+0.5,mean_profit_overstay_gap,'FaceColor',[0.2,0.3,0.3]); hold off;
    text(num_poles(end-1),mean_overstay_gap(end-1)-0.2,{'\leftarrow Lower is','     better'},'fontsize',13);
    text(num_poles(end),mean_profit_charging_gap(end)-0.2,{'\leftarrow Higher is','     better'},'fontsize',13);
    ylabel({'optimal-baseline ratio (%)'},'interpreter','latex');
    grid on; 
    title({'Sensitivity Analysis','with Number of Poles at Charging Station'});
    xlabel('# of poles','interpreter','latex'); 
    set(gca,'fontsize',15);
    legend([b2,b3,b1],{'profit_charging','profit_overstay','overstay'});
    
% style2: bar plot in each figure
elseif style == 2
    figure;
    h1 = subplot(211);
    b1 = bar(num_poles,mean_profit_gap);
    title({'Sensitivity Analysis','with Number of Poles at Charging Station'});
    set(gca,'fontsize',15); grid on;
    p1=get(h1,'position');
    legend('Profit')

    h2 = subplot(212,'position',[p1(1), p1(2)-p1(4)-0.1, p1(3), p1(4)]);
    b2 = bar(num_poles,-(mean_overstay_gap),'facecolor',[0.4660 0.6740 0.1880]); 
    xlabel('Number of poles'); 
    legend('Overstay');
    grid on; set(gca,'fontsize',15);
    p2=get(h2,'position');
    height=p1(2)+p1(4)-p2(2);
    h3=axes('position',[p2(1)-0.01 p2(2) p2(3) height],'visible','off');
    ylabel('optimal-baseline ratio (%)','visible','on','fontsize',18);
 
% style3: grouped bar chart 
elseif style == 3
    figure;
    b=bar(num_poles,[mean_profits_gap,-mean_overstay_gap]); % 09.20.2019 Teng: not change since not using
    b(1).FaceColor = [0 0.7 0.9];
    b(2).FaceColor = [0.4660 0.6740 0.1880];
    legend('Profit','Overstay');
    set(gca,'fontsize',15);
    ylabel('Optimal-baseline ratio (%)');
    xlabel('Number of poles');
    title({'Sensitivity Analysis','with Number of Poles at Charging Station'});
    grid on;
end

% saveas(gcf,'../recent-visualization/sensitivity_num_pole.png');
% saveas(gcf,'../recent-visualization/sensitivity_num_pole.fig');
% saveas(gcf,'../recent-visualization/sensitivity_num_pole','epsc');

%% visalize the quality of service
figure;
plot(num_poles,plot_service_tot,'linewidth',2); hold on;
plot(num_poles,plot_service_tot_base,'linewidth',2);
plot(num_poles,plot_service_tot_max,'-+b');plot(num_poles,plot_service_tot_min,'--b');
plot(num_poles,plot_service_tot_max_base,'-+r');plot(num_poles,plot_service_tot_min_base,'--r');
hold off; grid on;
% boxplot(box_service_tot,2:15); 
% boxplot(box_service_tot_base,2:15); hold off;grid on;
xlabel('Number of poles');
ylabel('Quality of service [0,1]');
set(gca, 'fontsize' , 15);
legend('W/ incentive','W/o incentive','location','nw'); 

%% visualize profit by charging option/overstay
figure(9);
subplot(121);
num_sim_in_monte = length(sim_results);
b=bar(num_poles, [profit_charging_c_sum, profit_charging_uc_sum, profit_overstay_sum]/num_sim_in_monte);
% xlabel('Number of poles');
ylabel('Profit ($)');
title({'Average Profit Distributions','[with Incentive Control]'});
xlabel('Number of poles');
legend('Charging-Flex','Charging-ASAP','Overstay');
% ylim([0, 200]);
set(gca,'fontsize',20); grid on;  
ylim([0 max(max([profit_charging_uc_sum, profit_charging_c_sum, profit_overstay_sum]/num_sim_in_monte))+10]);

subplot(122);
hb=bar(num_poles, [profit_charging_uc_sum_base, profit_overstay_sum_base]/num_sim_in_monte);
hb(1).FaceColor = [0.8500 0.3250 0.0980];
hb(2).FaceColor = [0.9290 0.6940 0.1250];
ylim([0 max(max([profit_charging_uc_sum, profit_charging_c_sum, profit_overstay_sum]/num_sim_in_monte))+10]);

xlabel('Number of poles');
% ylabel('Profit [$]');
title({'Average Profit Distributions','[without Incentive Control]'});
xlabel('Number of poles');
legend('Charging','Overstay')
% ylim([0, max()]);
set(gca,'fontsize',20); grid on; 
%
end