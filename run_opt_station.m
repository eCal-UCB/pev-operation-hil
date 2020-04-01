function opt = run_opt_station(station, k)
% optimization on station level
tic;
par = get_glob_par();
prb = get_glob_prb();

%% Functions
% fprintf('TEST: N_asap: %d\n',prb.N_asap);
% fprintf('TEST: length(TOU): %d\n',length(prb.TOU));

%% Run algorithm -- block coordinate descent

% compute time window of FLEX user
num_flex_user = length(station('FLEX_list'));
user_keys = station('FLEX_list');
flex_user_info = zeros(num_flex_user, 4);
for i = 1:num_flex_user
    user = station(user_keys{i,1});
    start_time = user.time.start/par.Ts; end_time = user.time.end/par.Ts;  
    N_remain = end_time-k/par.Ts;
    TOU_idx = k/par.Ts-start_time+1;
    flex_user_info(i, :) = [start_time, end_time, N_remain, TOU_idx]; % may consider to speed up by preallocate size
end
new_user = prb.user;
start_time = new_user.time;
existing_user_info = [start_time, -1, prb.N_flex, 1; flex_user_info]; % -1 indicate not yet specified, dimention (1 + num_flex_user) x 4

num_col_xk = max(existing_user_info(:, 3)); % set optimization time window
var_dim_constant = 2*num_col_xk+1;
xk = ones(var_dim_constant*(num_flex_user+1), 1); % [soc0_newuser .. socN_newuser u0_newuser .. uN_1_newuser; ...;  soc0_extuser .. socN_extsuser u0_extuser .. uN_1_extuser]; 
                                                  % - dimention var_dim_constant x (1+num_flex_user)

itermax = 1e4;
count = 0; improve = inf;
zk = ones(4,1);                         % [z_c, z_uc, y, 1];
vk = 1/3*ones(3,1);                     % [sm_c, sm_uc, sm_y];
Jk = zeros(itermax,1);
while count < itermax && improve >= 0 && abs(improve) >= par.opt.eps
    count = count + 1;
%     Jk(count) = J(zk,xk,vk);
    Jk(count) = constr_J(par,prb,zk,xk,vk,station,k,existing_user_info,var_dim_constant);
    
    % update init variables
    prb.z0 = zk; prb.x0 = xk; prb.v0 = vk; set_glob_prb(prb);
    
    % update control variables
    zk = argmin_z([],xk,vk);
    xk = argmin_x(zk,[],vk);
    vk = argmin_v(zk,xk,[]);
    
    % compute residual
    improve = Jk(count)-J(zk,xk,vk);
    
%     if mod(count,1) == 0
%         fprintf('[ OPT] iter: %d, improve: %.3f\n',count,improve);
%     end
end

opt.z = zk;
opt.tariff.flex = zk(1);
opt.tariff.asap = zk(2);
opt.tariff.overstay = zk(3);
opt.x = xk;
opt.flex.SOCs = xk(1:prb.N_flex+1);
opt.flex.powers = xk(prb.N_flex+2:end);
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