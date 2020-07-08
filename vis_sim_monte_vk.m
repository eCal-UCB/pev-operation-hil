function vis_sim_monte_vk(varargin)
% file_dir = 'monte-sim-results/';
% file_name = '06_03_20_21_14_monte_eps5_rand_seq_poles6.mat';
% load([file_dir file_name])
% load(fname);
if nargin == 0
    [fname, fpath] = uigetfile;
    data = load(fullfile(fpath,fname));
    sim_results = data.sim_results; 
%     sim_results_base = run_sim_baseline(data); 
elseif nargin == 1
    sim_results = varargin{1}; 
elseif nargin == 2
    sim_results = varargin{1}; 
    sim_results_base = varargin{2};
end
    
flex_prob = [];
asap_prob = [];
leave_prob = [];
for day = 1:length(sim_results)
    sim = sim_results{day};
    choice_probs = sim.choice_probs(~isnan(sim.choice_probs(:,1)),:);
    flex_prob = [flex_prob mean(choice_probs(:,1))];
    asap_prob = [asap_prob mean(choice_probs(:,2))];
    leave_prob = [leave_prob mean(choice_probs(:,3))];
end

figure;
plot(flex_prob, 'LineWidth', 2)
hold on 
plot(flex_prob + asap_prob, 'LineWidth', 2)
hold off
xlabel('Number of simulations', 'FontSize', 16)
ylabel('Probability', 'FontSize', 16)
legend('AVG. FLEX Prob', 'AVG. Prob (Cum. of FLEX and ASAP)', 'FontSize', 12, 'Location', 'Best')
title('Average Choice Probabilities (Controlled)', 'FontSize', 16)

% flex_prob = [];
% asap_prob = [];
% leave_prob = [];
% for day = 1:length(sim_results_base)
%     sim = sim_results_base{day};
%     choice_probs = sim.choice_probs(~isnan(sim.choice_probs(:,1)),:);
%     flex_prob = [flex_prob mean(choice_probs(:,1))];
%     asap_prob = [asap_prob mean(choice_probs(:,2))];
%     leave_prob = [leave_prob mean(choice_probs(:,3))];
% end
% 
% figure;
% plot(flex_prob, 'LineWidth', 2)
% hold on 
% plot(flex_prob + asap_prob, 'LineWidth', 2)
% hold off
% xlabel('Number of simulations')
% ylabel('Probability')
% legend('AVG. FLEX Prob', 'AVG. Prob (Cum. of FLEX and ASAP)', 'FontSize', 12, 'Location', 'Best')
% title('Average Choice Probabilities (baseline)', 'FontSize', 16)
end