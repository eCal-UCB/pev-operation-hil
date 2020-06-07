function [station, opt] = run_opt_station(station, k)
% optimization on station level
tic;
par = get_glob_par();
prb = get_glob_prb();

%% Functions
% fprintf('TEST: N_asap: %d\n',prb.N_asap);
% fprintf('TEST: length(TOU): %d\n',length(prb.TOU));

%% Run algorithm -- block coordinate descent

test_tic = tic;
for idx = 1:1
%     disp(idx)
% compute time window of FLEX user
num_flex_user = length(station('FLEX_list'));
user_keys = station('FLEX_list');
flex_user_info = zeros(num_flex_user, 6);
for i = 1:num_flex_user
    user = station(user_keys{i});
    start_time = user.time.start/par.Ts; end_time = user.time.end/par.Ts;  
    N_remain = end_time-k/par.Ts;
    TOU_idx = k/par.Ts-start_time+1; % indicate electricity price, also for calculating updated SOC
    SOC_need = user.prb.user.SOC_need;
    try
        SOC_now = user.prb.user.SOC_init + sum(user.powers(1:TOU_idx-1)) * user.par.eff * par.Ts / user.prb.user.batt_cap;
    catch
        SOC_now = user.prb.user.SOC_init + sum(user.powers(1:end)) * user.par.eff * par.Ts / user.prb.user.batt_cap;
    end
    flex_user_info(i, :) = [start_time, end_time, N_remain, TOU_idx, SOC_need, SOC_now];
end
new_user = prb.user;
start_time = new_user.time/par.Ts;
existing_user_info = [start_time, -1, prb.N_flex, 1, new_user.SOC_need, new_user.SOC_init; flex_user_info]; % -1 indicate not yet specified, dimention (1 new user + num_flex_user) x 4

num_col_xk = max(existing_user_info(:, 3)); % set optimization time window, T_end^FLEX
var_dim_constant = 2*num_col_xk+1;
xk = ones(var_dim_constant*(num_flex_user+2), 1); % [soc0_newuser .. socN_newuser u0_newuser .. uN_1_newuser; ...;  soc0_extuser .. socN_extsuser u0_extuser .. uN_1_extuser]; 
                                                  % - dimention
                                                  % var_dim_constant x
                                                  % (1+num_flex_user) --->
                                                  % reshape to one column
                                                  % (2+num_flex_user) for
                                                  % demand charge

itermax = 1e4;
count = 0; improve = inf;
zk = ones(4,1);                         % [z_c, z_uc, y, 1];
vk = [0.45 0.45 0.1]';                     % [sm_c, sm_uc, sm_y];
Jk = zeros(itermax,1);
while count < itermax && improve >= 0 && abs(improve) >= par.opt.eps
    count = count + 1;
%     Jk(count) = J(zk,xk,vk);
    Jk(count) = constr_J(par,prb,zk,xk,vk,station,k,existing_user_info,var_dim_constant);
    
    % update init variables
    prb.z0 = zk; prb.x0 = xk; prb.v0 = vk; set_glob_prb(prb);
    % update control variables - station
    zk = argmin_z_station([],xk,vk,station,k,existing_user_info,var_dim_constant);
    xk = argmin_x_station(zk,[],vk,station,k,existing_user_info,var_dim_constant);
    vk = argmin_v_station(zk,xk,[],station,k,existing_user_info,var_dim_constant);

    % compute residual
    improve = Jk(count)-constr_J(par,prb,zk,xk,vk,station,k,existing_user_info,var_dim_constant);
    
    if count > 0
%         fprintf('[ OPT] iter: %d, improve: %.3f, objective: %.3f\n',count,improve,Jk(count));
        fprintf('[ OPT] iter: %d, probability: ', count);
        vk
    end
    
%     if mod(count,1) == 0
%         fprintf('[ OPT] iter: %d, improve: %.3f\n',count,improve);
%     end
end
end
fprintf('elapsed time for running %.f times is: %.2f seconds. \n', idx, toc(test_tic))

% ===== TODO: update station users profile
% iterate thru EV in station
% replace x and powers in user
% add field of SOC in user

% update flex user profile
for i = 1:length(user_keys)
    user = station(user_keys{i});
    end_time = user.time.end/par.Ts;  
    N_remain = end_time-k/par.Ts;
    user.x = xk((i)*var_dim_constant+1:(i+1)*var_dim_constant);
    user.powers = xk((i)*var_dim_constant+N_remain+2:(i+1)*var_dim_constant);
    user.SOC = xk((i)*var_dim_constant+1:i*var_dim_constant+N_remain+1);
    station(user_keys{i}) = user;
end

zk

opt.z = zk;
opt.tariff.flex = zk(1);
opt.tariff.asap = zk(2);
opt.tariff.overstay = zk(3);
opt.x = xk; % for all flex user/vehicle
opt.peak_pow = xk(end);
opt.flex.SOCs = xk(1:prb.N_flex+1); % record new user flex
opt.flex.powers = xk(num_col_xk+2:var_dim_constant); % record new user asap
opt.asap.powers = prb.station.pow_max;
opt.v = vk;
opt.prob.flex = vk(1);
opt.prob.asap = vk(2);
opt.prob.leave = vk(3);
opt.J = Jk(1:count);
opt.num_iter = count;
opt.prb = prb;
opt.par = par;

opt.time.start = prb.user.time;
opt.time.end_flex = prb.user.time + prb.user.duration;
opt.time.end_asap = prb.user.time + prb.N_asap*par.Ts;

fprintf('[%s OPT] DONE (%.2f sec) sum(vk) = %.2f, iterations = %d\n',datetime('now'),toc,sum(vk),count);
end