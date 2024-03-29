function xk = argmin_x_station(z,~,v,station,k,existing_user_info,var_dim_constant)
% flex_user_info(i, :) = [start_time, end_time, N_remain, TOU_idx, SOC_need, SOC_now];

par = get_glob_par();
prb = get_glob_prb();

if (v < 0) | (sum(v) < 1-par.soft_v_eta) | (sum(v) > 1+par.soft_v_eta)
    warning('[ Warning] invalid $v$');
end

% cost function
J = @(xk) constr_J(xk);
num_flex_user = size(existing_user_info,1);

% constraints
A_ineq = []; b_ineq = zeros(num_flex_user,1); % inequality constraint
lb = zeros(var_dim_constant*num_flex_user, 1); % lower bound - power min
ub = zeros(var_dim_constant*num_flex_user, 1); % upper bound - power max
A_eq = []; b_eq = []; % equality constraints - system dynamics

% for all existing + new users
for idx = 1:size(existing_user_info,1)
   % ==========================
   % get params
   % ==========================
   N = existing_user_info(idx,3); % get user N_flex, aka the remain duration
   N_max = (var_dim_constant-1)/2;
   SOC_feas_max = (prb.station.pow_max * N * par.Ts * par.eff)/prb.user.batt_cap;
   SOC_need = existing_user_info(idx,5);
   SOC_init = existing_user_info(idx,6);
   SOC_need = min(SOC_need, SOC_feas_max+SOC_init);
   
   % ======================================================================
   % inequality constraint -- charging requirements, power bounds
   % ======================================================================
   % SOC requirement
   A_ineq_user = [zeros(1,N) -1 zeros(1,var_dim_constant-N-1)];
   A_ineq = blkdiag(A_ineq, A_ineq_user);
   b_ineq(idx) = -SOC_need; % SOC_need
   
   % lower bound - energy/power min
   lb1 = [zeros(N,1); SOC_need; zeros(N_max-N,1)]; % <NOTE> this is redundant with inequality constraint
   lb2 = [prb.station.pow_min.*ones(N,1); zeros(N_max-N,1)];
   lb((idx-1)*var_dim_constant+1:idx*var_dim_constant,1) = [lb1;lb2]; %1:soc, 2:power
   
   % upper bound - energy/power max
   ub1 = [ones(N,1); 1; zeros(N_max-N,1)]; 
   ub2 = [prb.station.pow_max.*ones(N,1); zeros(N_max-N,1)];
   ub((idx-1)*var_dim_constant+1:idx*var_dim_constant,1) = [ub1;ub2];
   
   % ======================================
   % equality constraint -- system dynamics
   % ======================================
   % initial soc
   C1L = [1 zeros(1,N)];
   C1R = zeros(1,var_dim_constant-N-1); 
   d1 = SOC_init; 
   
   % system dynamics
   C2L = [diag(-1.*ones(1,N)) zeros(N,1) zeros(N,N_max-N)] + [zeros(N,1)  diag(ones(1,N)) zeros(N,N_max-N)];
   C2R = [-diag(par.eff*par.Ts/prb.user.batt_cap*ones(1,N)) zeros(N,N_max-N)];
   d2 = zeros(N,1);
   

   A_eq = blkdiag(A_eq, [C1L C1R; C2L C2R]);
   b_eq = [b_eq; [d1; d2]];
end

% solve optimization
options = optimoptions('fmincon','Display','off');
% options.Algorithm = 'sqp';

xk = fmincon(J,ones(size(A_ineq,2),1),A_ineq,b_ineq,A_eq,b_eq,lb,ub,[],options);

function J = constr_J(x)
    N_max = (var_dim_constant-1)/2;
% par, prb, z, x, v
% station - container.maps object 
% k = global time index
    % existing flex user
    user_keys = station('FLEX_list');
    existing_flex_obj = 0;
    for i = 2:size(existing_user_info,1) % sum of users
        adj_constant = (i-1) * var_dim_constant; % constant to identify where to start on x
        duration = existing_user_info(i,3); TOU_idx = existing_user_info(i,4);
        user = station(user_keys{1,i-1});
        overstay_cost = (user.time.leave - user.time.end) * user.z(3);
%         existing_flex_obj = existing_flex_obj + (sum(x(adj_constant+duration+2:adj_constant+2*duration+1,1).*(user.prb.TOU(TOU_idx:end) - user.price)) - overstay_cost);
        existing_flex_obj = existing_flex_obj + (sum(x(adj_constant+N_max+2:adj_constant+N_max+2+duration-1,1).*(user.prb.TOU(TOU_idx:end) - user.price)) - overstay_cost);
    end
    % existing asap user
    user_keys = station('ASAP_list');
    existing_asap_obj = 0;
    for i = 1:length(user_keys) % sum of users
        user = station(user_keys{1,i});
        overstay_cost = (user.time.leave - user.time.end) * user.z(3);
        TOU_idx = (k-user.time.start)/par.Ts+1;
        existing_asap_obj = existing_asap_obj + (sum(mean(user.asap.powers)*(user.prb.TOU(TOU_idx:end) - user.price)) - overstay_cost);
    end
    
    %%% planned asap power profile %%%
    asap_power_sum_profile = zeros(1,var_dim_constant);
    it = 0;
    for t = k : k + (N_max-1)*par.Ts
        it = it + 1;
        for i = 1:length(station('ASAP_list'))
            opt = station(user_keys{1,i});
            if k <= opt.time.end
                asap_power_sum_profile(it) = asap_power_sum_profile(it) + interp1(opt.time.start:par.Ts:opt.time.end-par.Ts,opt.powers,k);
            end
        end
    end
    
    % ==== missing demand charge ====
    
    % part 1: case 1 - charging-FLEX
%     new_flex_obj = (sum((x(prb.N_flex+2:2*prb.N_flex+1).*(prb.TOU(1:prb.N_flex) - z(1)))...
%                         +par.lambda.x .* x(prb.N_flex+2:2*prb.N_flex+1))...
%                     +par.lambda.z_c*z(1)^2)...
%                   +par.lambda.h_c * 1/z(3); % with convergence regularization
              
%     new_flex_obj = sum((x(prb.N_flex+2:2*prb.N_flex+1).*(prb.TOU(1:prb.N_flex) - z(1)))...
%                     +par.lambda.z_c*z(1)^2)...
%                   +par.lambda.h_c * 1/z(3); % with convergence regularization
    new_flex_obj = sum((x(N_max+2:N_max+2+prb.N_flex-1).*(prb.TOU(1:prb.N_flex) - z(1)))...
                    +par.lambda.z_c*z(1)^2)...
                  +par.lambda.h_c * 1/z(3); % with convergence regularization

%     new_flex_obj = (sum((x(prb.N_flex+2:2*prb.N_flex+1,1).*(prb.TOU(1:prb.N_flex) - z(1)))...
%                         +par.lambda.x .* x(prb.N_flex+2:2*prb.N_flex+1,1)))...
%                   +par.lambda.h_c * 1/z(3); % without convergence regularization 
              
    % part 2: case 2 - charging-ASAP
    new_asap_obj = (sum((prb.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2))))...
                    +par.lambda.z_c*z(2)^2)...
                  +par.lambda.h_c * 1/z(3); % with convergence regularization

%     new_asap_obj = sum(prb.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2)))...
%                   +par.lambda.h_uc * 1/z(3); % without convergence regularization
%               +par.lambda.x .* x(prb.N_asap+2:2*prb.N_asap+1,1)
              
    % part 3: case 3 - leave
    new_leave_obj = sum(prb.station.pow_max*(prb.TOU(1:prb.N_asap) - 0));
    
%     J = dot([new_flex_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * x(end); 
%         new_asap_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * x(end); 
%         new_leave_obj], v);
    J = dot([new_flex_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * max(asap_power_sum_profile + sum(reshape(x,var_dim_constant,length(x)/var_dim_constant)',1));
        new_asap_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * max(asap_power_sum_profile + [ones(1,prb.N_asap)*prb.station.pow_max zeros(1,var_dim_constant-prb.N_asap)]); 
        new_leave_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * max(asap_power_sum_profile)], v);
%     J = dot([new_flex_obj; 
%         new_asap_obj; 
%         new_leave_obj], v);
end
end