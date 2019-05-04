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
        figure; num_subplot = 4; count = 1;
        subplot(eval([num2str(num_subplot) '1' num2str(count)])); count = count + 1;
        plot(sim.t,sim.power,'linewidth',1.5); xlim([sim.t(1) sim.t(end)]);
        grid on; xlabel('hour of the day'); ylabel('power (kW)');
        set(gca,'fontsize',15);

        subplot(eval([num2str(num_subplot) '1' num2str(count)])); count = count + 1;
        plot(sim.t,sim.profit,'linewidth',1.5); hold on;
        plot(sim.t,cumsum(sim.profit),'linewidth',1.5); hold off; xlim([sim.t(1) sim.t(end)]);
        grid on; xlabel('hour of the day'); ylabel('profit ($)');
        legend('instant','net')
        set(gca,'fontsize',15);

        subplot(eval([num2str(num_subplot) '1' num2str(count)])); count = count + 1;
        plot(sim.t,sim.occ,'linewidth',1.5); hold on;
        plot(sim.t,sim.overstay,'linewidth',1.5);
        plot(sim.t,sim.charging,'linewidth',1.5); hold off;
        xlim([sim.t(1) sim.t(end)]); legend('total','overstay','charging');
        grid on; xlabel('hour of the day'); ylabel('# of vehicles');
        set(gca,'fontsize',15); 

        subplot(eval([num2str(num_subplot) '1' num2str(count)])); count = count + 1;
        plot(sim.t,cumsum(sim.overstay_duration),'linewidth',1.5); xlim([sim.t(1) sim.t(end)]);
        grid on; set(gca,'fontsize',15); xlabel('hour of the day'); ylabel([{'overstay duration'}, {'(hours)'}]);
    end
    
    if options.choices
        figure; 
        choices = sim.choice(~isnan(sim.choice));
        choice_probs = sim.choice_probs(sim.choice_probs(:,1)~=0,:);
        choice_labels = zeros(length(choices),1);
        for j = 1:length(choice_labels)
            choice_labels(j) = (choices(j)==0)*1/2*choice_probs(j,1) ...
                               + (choices(j)==1)*(choice_probs(j,1)+1/2*choice_probs(j,2))...
                               + (choices(j)==2)*(choice_probs(j,1)+choice_probs(j,2)+1/2*choice_probs(j,3));
        end
        area(1:length(find(sim.choice_probs(:,1))), sim.choice_probs(sim.choice_probs(:,1)~=0,:)); hold on;
        text(1.5,choice_probs(1,1)-0.04,'charging flex','fontsize',15,'color','white');
        text(1.5,sum(choice_probs(1,1:2))-0.04,'charging asap','fontsize',15,'color','white');
        text(1.5,sum(choice_probs(1,1:3))-0.04,'leaving without charging','fontsize',15,'color','white');
        h = scatter(1:length(choices), choice_labels,'k','filled'); hold off;
        xlabel('event#'); ylabel('probability [0,1]'); 
        title([{'Choice Probabilities and Decisions '},...
               {sprintf('(total: %d, flex: %d, asap: %d, leave: %d)',...
                        length(choices),...
                        sum(choices==0),...
                        sum(choices==1),...
                        sum(choices==2)...
                        )}]);
        xlim([1,length(find(sim.choice_probs(:,1)))]); ylim([0,1]);
        legend(h,{'choice'})
        set(gca,'fontsize',15);
    end
end
end