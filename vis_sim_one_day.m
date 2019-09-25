function options = vis_sim_one_day(varargin)
% options include: 
if nargin == 0 && nargout == 1
    options.display = true;
    options.temporals = true;
    options.choices = true;
    return;
end

if nargout == 0 % visualize
    if nargin == 1 
        options.display = true;
        options.temporals = true;
        options.choices = true;
    elseif nargin == 2
        options.display = varargin{2}.display;
        options.temporals = varargin{2}.temporals;
        options.choices = varargin{2}.choices;
    end
end

sim = varargin{1};
if options.display
    if options.temporals
        figure('position',[1,1,640,704]); num_subplot = 5; count = 1;
        
        % power consumption
        subplot(eval([num2str(num_subplot) '1' num2str(count)])); count = count + 1;
        plot(sim.t,sim.power,'linewidth',1.5); xlim([sim.t(1) sim.t(end)]);
        grid on; xlabel('hour of the day'); ylabel('power (kW)');
        set(gca,'fontsize',15);
        
        % profit
        subplot(eval([num2str(num_subplot) '1' num2str(count)])); count = count + 1;
%         plot(sim.t,sim.profit_charging_uc,'linewidth',1.5); hold on;
%         plot(sim.t,sim.profit_charging_c,'linewidth',1.5); hold on;
%         plot(sim.t,sim.profit_overstay,'linewidth',1.5); hold on;
%         plot(sim.t,cumsum(sim.profit_charging_uc),'linewidth',1.5); hold on;
%         plot(sim.t,cumsum(sim.profit_charging_c),'linewidth',1.5); hold on;
%         plot(sim.t,cumsum(sim.profit_overstay),'linewidth',1.5); hold off; xlim([sim.t(1) sim.t(end)]);
%         grid on; xlabel('hour of the day'); ylabel('profit ($)');
%         legend('instant-charging-uc', 'instant-charging-c', 'instant-overstay', ...
%             'net-charging-uc', 'net-charging-c', 'net-overstay')
        plot(sim.t,sim.profit_charging_uc+sim.profit_charging_c+sim.profit_overstay,'linewidth',1.5); hold on;
        plot(sim.t,cumsum(sim.profit_charging_uc+sim.profit_charging_c+sim.profit_overstay),'linewidth',1.5); 
        hold off; xlim([sim.t(1) sim.t(end)]);
        grid on; xlabel('hour of the day'); ylabel('profit ($)');
        legend('Instant','Net');
        set(gca,'fontsize',15);

        % occupancy
        subplot(eval([num2str(num_subplot) '1' num2str(count)])); count = count + 1;
        plot(sim.t,sim.occ.total,'linewidth',1.5); hold on;
        plot(sim.t,sim.occ.overstay,'linewidth',1.5);
        plot(sim.t,sim.occ.charging,'linewidth',1.5); hold off;
        xlim([sim.t(1) sim.t(end)]); legend('total','overstay','charging');
        grid on; xlabel('hour of the day'); ylabel('# of vehicles');
        set(gca,'fontsize',15); 

        % overstay duration
        subplot(eval([num2str(num_subplot) '1' num2str(count)])); count = count + 1;
        plot(sim.t,cumsum(sim.overstay_duration),'linewidth',1.5); xlim([sim.t(1) sim.t(end)]);
        grid on; set(gca,'fontsize',15); xlabel('hour of the day'); ylabel([{'overstay duration'}, {'(hours)'}]);
        
        % number of serviced vehicles
        subplot(eval([num2str(num_subplot) '1' num2str(count)])); count = count + 1;
        stem(sim.t,sim.num_service,'linewidth',1.5,'markersize',eps); hold on;
        plot(sim.t,cumsum(sim.num_service),'linewidth',1.5); hold off;
        xlim([sim.t(1) sim.t(end)]); legend('instant','net');
        grid on; set(gca,'fontsize',15); xlabel('hour of the day'); ylabel([{'# of service'}]);
        
        
        
        saveas(gcf,'recent-visualization/one_day_temporal.png');
        saveas(gcf,'recent-visualization/one_day_temporal.fig');
        saveas(gcf,'recent-visualization/one_day_temporal','epsc');
    end
    
    if options.choices
        figure; 
        choices = sim.choice(~isnan(sim.choice));
        choice_probs = sim.choice_probs(~isnan(sim.choice_probs(:,1)),:);
        choice_labels = zeros(length(choices),1);
        choice_times = sim.events.time(logical(sim.events.triggered));
        choice_times_str = num2str(round(choice_times*100)/100);
        choice_times_str_filtered = cell(length(choice_labels),1);
        for i = 1:length(choice_labels)
            choice_times_str_filtered{i} = ['(' strtrim(choice_times_str(i,:)) ')'];
        end
        
        for j = 1:length(choice_labels)
            choice_labels(j) = (choices(j)==0)*1/2*choice_probs(j,1) ...
                               + (choices(j)==1)*(choice_probs(j,1)+1/2*choice_probs(j,2))...
                               + (choices(j)==2)*(choice_probs(j,1)+choice_probs(j,2)+1/2*choice_probs(j,3));
        end
        area(1:length(choice_probs), choice_probs); hold on;
%         text(1.5,choice_probs(1,1)-0.04,'charging flex','fontsize',15,'color','white');
        text(1.5,0.04,'charging flex','fontsize',15,'color','white');
        text(1.5,sum(choice_probs(1,1:2))-0.04,'charging asap','fontsize',15,'color','white');
        text(1.5,sum(choice_probs(1,1:3))-0.04,'leaving without charging','fontsize',15,'color','white');
        
        h = scatter(1:length(choices), choice_labels,'k','filled'); hold off;
        text(1:length(choice_probs), choice_labels+0.01, choice_times_str_filtered, 'fontsize',13);
        
        xlabel('event#'); ylabel('probability [0,1]'); 
        title([{'Choice Probabilities and Decisions '},...
               {sprintf('(total: %d, flex: %d, asap: %d, leave: %d)',...
                        length(choices),...
                        sum(choices==0),...
                        sum(choices==1),...
                        sum(choices==2)...
                        )}]);
        xlim([1,length(choice_probs)]); ylim([0,1]);
        legend(h,{'choice (hours)'})
        set(gca,'fontsize',15);
        %%
        saveas(gcf,'recent-visualization/one_day_choice.png');
        saveas(gcf,'recent-visualization/one_day_choice.fig');
        saveas(gcf,'recent-visualization/one_day_choice','epsc');
    end
end
end