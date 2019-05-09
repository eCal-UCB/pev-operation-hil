function vis_sensitivity_analysis(varargin)
if nargin == 0
    [fname, fpath] = uigetfile;
    data = load(fullfile(fpath,fname));
    monte_results = data.monte_results; 
elseif nargin == 1
    monte_results = varargin{1};
end

% Data process
num_monte = length(monte_results);

num_poles = monte_results{1}.optimal{1}.par.sens_analysis.num_poles;
mean_profits_gap = zeros(num_monte,1);
mean_overstay_gap = zeros(num_monte,1);
for i = 1:length(monte_results)
    sim_results = monte_results{i}.optimal;
    sim_results_base = monte_results{i}.baseline;
    
    num_sims = sim_results{1}.par.monte.num_sims;
    
    % with controller
    overstay_tot = zeros(num_sims,1);
    overstay_mean = zeros(num_sims,1);
    overstay_penalty_mean = zeros(num_sims,1);
    profit = zeros(num_sims,1);
    service_tot = zeros(num_sims,1);
    for n = 1:num_sims
        overstay_mean(n) = mean(sim_results{n}.overstay_duration(sim_results{n}.overstay_duration~=0));
        overstay_penalty_mean(n) = mean(sim_results{n}.control(sim_results{n}.control(:,3)~=0,3));
        overstay_tot(n) = sum(sim_results{n}.overstay_duration);
        profit(n) = sum(sim_results{n}.profit);
        service_tot(n) = sum(sim_results{n}.num_service);
    end

    % without controller
    overstay_tot_base = zeros(num_sims,1);
    overstay_mean_base = zeros(num_sims,1);
    overstay_penalty_mean_base = zeros(num_sims,1);
    profit_base = zeros(num_sims,1);
    service_tot_base = zeros(num_sims,1);
    for n = 1:num_sims
        overstay_mean_base(n) = mean(sim_results_base{n}.overstay_duration(sim_results_base{n}.overstay_duration~=0));
        overstay_penalty_mean_base(n) = mean(sim_results_base{n}.control(sim_results_base{n}.control(:,3)~=0,3));
        overstay_tot_base(n) = sum(sim_results_base{n}.overstay_duration);
        profit_base(n) = sum(sim_results_base{n}.profit);
        service_tot_base(n) = sum(sim_results_base{n}.num_service);
    end
    
    % compute the improvements
    mean_profits_gap(i) = (mean(profit)/mean(profit_base)-1)*100;
    mean_overstay_gap(i) = (mean(overstay_mean)/mean(overstay_mean_base)-1)*100;
end


%% Visualization
close;figure; 
% scatter(num_poles,mean_profits_gap,100*ones(size(num_poles)),'filled'); hold on;
% scatter(num_poles,mean_overstay_gap,100*ones(size(num_poles)),'filled'); hold off;
 
b1=bar(num_poles,mean_overstay_gap,'facecolor',[0.4660 0.6740 0.1880]); hold on;
b2=bar(num_poles,mean_profits_gap,'FaceColor',[0,0.7,0.9]); hold off;
text(num_poles(end-1),mean_overstay_gap(end-1)-0.2,{'\leftarrow Lower is','     better'},'fontsize',13);
text(num_poles(end),mean_profits_gap(end)-0.2,{'\leftarrow Higher is','     better'},'fontsize',13);
ylabel({'optimal-baseline ratio (%)'},'interpreter','latex');
grid on; 
title({'Sensitivity Analysis','with Number of Poles at Charging Station'});
xlabel('# of poles','interpreter','latex'); 
set(gca,'fontsize',15);
legend([b2,b1],{'profit','overstay'});

saveas(gcf,'recent-visualization/sensitivity_num_pole.png');
saveas(gcf,'recent-visualization/sensitivity_num_pole.fig');
saveas(gcf,'recent-visualization/sensitivity_num_pole','epsc');
end