% This script is to simulate EV charging station operations where the
% charging tariff is determined real-time with taking account into EV
% drivers' behaviors. The overall objective of the tariff control is to
% minimize the operation cost of the system operater. 
%
% There are three choices that each EV driver can make at arrival:
% (i) charging with flexibility, (ii) charging as soon as possible, and
% (iii) leaving without charging.
%
% At each arrival of EV driver, user specific parameters, e.g., battery
% capacity, desired parking durations, initial SOC, and needed SOC level
% (for next mobility demand) are randomly sampled from an empirical
% probability distribution function that is generated with TELD dataset 
% [X].
%
% THIS WORK IS A PART OF EE227C COURSE PROJECT AT UC BERKELEY.
% last modified, May 2019.
%
% Contributors: Sangjae Bae, Teng Zeng, Bertrand Travacca.

clear
%% Initialization
disp('[ INIT] initializing...');
par = set_glob_par(init_params());
disp('[ INIT] DONE');


%% Simulation
station = containers.Map; 
station('num_occupied_pole') = 0;
tot_visit = 0;
tot_power = zeros(size(par.sim.starttime:par.Ts:par.sim.endtime));
tot_rev  = zeros(size(par.sim.starttime:par.Ts:par.sim.endtime));
tot_occ  = zeros(size(par.sim.starttime:par.Ts:par.sim.endtime));
i = 0;
for k = par.sim.starttime:par.Ts:par.sim.endtime
    i = i+1;
    
    % random visit
    rv = rand;
    if rv <= interp1(0:23,par.pdf.visit,k)
        tot_visit = tot_visit + 1;
       % TODO: sample values from pdf
       inp.time = k;
       inp.SOC_init = 0.2;
       inp.SOC_need = 0.4; % add infeasible scenario
       inp.batt_cap = 24;
       inp.duration = 3;
       set_glob_prb(init_prb(inp));
       
       if inp.duration >= par.sim.endtime - k
           break;
       end
       
       % find optimal tariff
       opt = run_opt();
       
       % driver makes choice
       rc = rand;
       if rc <= opt.prob.flex
           opt.choice = 0; % charging with flexibility
           opt.time.end = opt.time.end_flex;
           opt.powers = opt.flex.powers;
           opt.price = opt.flex.tariff;
       elseif rc <= opt.prob.flex + opt.prob.asap
           opt.choice = 1; % charging as soon as possible
           opt.time.end = opt.time.end_asap;
           opt.powers = opt.asap.powers;
           opt.price = opt.asap.tariff;
       else
           opt.choice = 2; % leaving without charging
       end
       fprintf('[ EVENT] time = %.2f, CHOICE = %s\n',k,par.dcm.choices{opt.choice+1});
       if opt.choice == 0 || opt.choice == 1
           station('num_occupied_pole') = station('num_occupied_pole') + 1;
           station(['EV' num2str(tot_visit)]) = opt;
       end
    end
    
    % update power sum
    keys = station.keys();
    if ~isempty(keys)
        for ev = keys
            if contains(ev{1},'EV')
                if  k <= station(ev{1}).time.end
                    TOU = interp1(0:23,par.TOU,k,'nearest');
                    if length(station(ev{1}).powers) > 1
                        power = interp1(linspace(station(ev{1}).time.start, ...
                                            station(ev{1}).time.end,...
                                            length(station(ev{1}).powers)), ...
                                            station(ev{1}).powers, k);
                    elseif length(station(ev{1}).powers) == 1
                        power = station(ev{1}).powers;
                    end
                    
                    tot_power(i) = tot_power(i) + power;
                    tot_rev(i) = tot_rev(i) + par.Ts * power * (station(ev{1}).price - TOU);
                else
                    % TB Removed: assume the vehicle does not overstay
                    
                    % TODO: overstay scenario
                    % - check overstaying (generate random var)
                    % - if stay, charge fee
                    % - if not stay, remove the key
                    station.remove(ev{1});
                    station('num_occupied_pole') = station('num_occupied_pole') - 1;
                end
            elseif contains(ev{1},'occ')
                tot_occ(i) = station('num_occupied_pole');
            end
        end
    end
end


%% Visualization
t = par.sim.starttime:par.Ts:par.sim.endtime;
figure(1); 
subplot(311);
plot(t,tot_power,'linewidth',1.5);
grid on; xlabel('hour of the day'); ylabel('kW'); title('Power Consumption');
set(gca,'fontsize',15);
subplot(312);
plot(t,tot_rev,'linewidth',1.5); hold on;
plot(t,cumsum(tot_rev),'linewidth',1.5); hold off;
grid on; xlabel('hour of the day'); ylabel('$'); title('Revenue');
set(gca,'fontsize',15);
subplot(313);
plot(t,tot_occ,'linewidth',1.5);
grid on; xlabel('hour of the day'); ylabel('#'); title('Occupancy');
set(gca,'fontsize',15);
