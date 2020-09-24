function zk = argmin_z_station(~,x,v,station,k,existing_user_info,var_dim_constant)
par = get_glob_par();
prb = get_glob_prb();

if (v < 0) | (sum(v) < 1-par.soft_v_eta) | (sum(v) > 1+par.soft_v_eta)
    error('[ ERROR] invalid $v$');
end

% LSE function
lse = @(z) log(exp(dot(prb.THETA(1,:),z))+exp(dot(prb.THETA(2,:),z))+exp(dot(prb.THETA(3,:),z)));

% cost function
% J = @(z) dot([sum((x(prb.N_flex+2:end).*(prb.TOU(1:prb.N_flex) - z(1))).^2) + par.lambda.h_c * 1/z(3);
%             sum((par.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2))).^2) + par.lambda.h_uc * 1/z(3);
%             1/3*sum((par.station.pow_max*(prb.TOU(1:prb.N_asap) - 0)).^2)],v) ...
%          + par.mu * (lse(z) - z' * prb.THETA' * v);
% J = @(z) dot([(sum((x(prb.N_flex+2:end).*(prb.TOU(1:prb.N_flex) - z(1)))+par.lambda.x.*x(prb.N_flex+2:end))+par.lambda.z_c*z(1)^2) + par.lambda.h_c * 1/z(3);
%             sum((prb.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2)))+par.lambda.z_uc*z(2)^2) + par.lambda.h_uc * 1/z(3);
%             sum(prb.station.pow_max*(prb.TOU(1:prb.N_asap) - 0))],v) ...
%          + par.mu * (lse(z) - z' * prb.THETA' * v);

J = @(zk) constr_J(zk);

% inequality constratins
A = [1 -1 0 0]; b = 0; % charging tariff for charging asap must be bigger
     
% lower and upper bound
% lb = [min(prb.TOU) * ones(2,1); 2*min(prb.TOU); 0];
% ub = [2 * max(prb.TOU) * ones(2,1); 10 * max(par.TOU); 1];
lb = [max(prb.TOU) * ones(3,1); 0];
ub = [2 * max(prb.TOU) * ones(2,1); 10 * max(prb.TOU); 1];

% equality constraint
Aeq = [0 0 0 1]; beq = 1;

% solve optimization
options = optimoptions('fmincon','Display','off');
% options.Algorithm = 'sqp';
zk = fmincon(J,prb.z0,A,b,Aeq,beq,lb,ub,[],options);

function J = constr_J(z)
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
        existing_flex_obj = existing_flex_obj + (sum(x(adj_constant+duration+2:adj_constant+2*duration+1,1).*(user.prb.TOU(TOU_idx:end) - user.price)) - overstay_cost);
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
    N_max = (var_dim_constant-1)/2;
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
    
    % part 1: case 1 - charging-FLEX
    new_flex_obj = (sum((x(prb.N_flex+2:2*prb.N_flex+1,1).*(prb.TOU(1:prb.N_flex) - z(1)))...
                        +par.lambda.x .* x(prb.N_flex+2:2*prb.N_flex+1,1))...
                    +par.lambda.z_c*z(1)^2)...
                  +par.lambda.h_c * 1/z(3); % with convergence regularization

              
              
%     new_flex_obj = (sum((x(prb.N_flex+2:2*prb.N_flex+1,1).*(prb.TOU(1:prb.N_flex) - z(1)))...
%                         +par.lambda.x .* x(prb.N_flex+2:2*prb.N_flex+1,1)))...
%                   +par.lambda.h_c * 1/z(3); % without convergence regularization 
              
    % part 2: case 2 - charging-ASAP
    new_asap_obj = (sum((prb.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2)))...
                    +par.lambda.z_c*z(2)^2))...
                  +par.lambda.h_c * 1/z(3); % with convergence regularization

%     new_asap_obj = (sum((prb.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2)))...
%                         +par.lambda.x .* x(prb.N_asap+2:2*prb.N_asap+1,1)))...
%                   +par.lambda.h_uc * 1/z(3); % without convergence regularization
              
    % part 3: case 3 - leave
    new_leave_obj = sum(prb.station.pow_max*(prb.TOU(1:prb.N_asap) - 0));
    
%     J = dot([new_flex_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * x(end); 
%         new_asap_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * x(end); 
%         new_leave_obj], v) + par.mu * (lse(z) - z' * prb.THETA' * v);
    J = dot([new_flex_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * max(asap_power_sum_profile + sum(reshape(x,var_dim_constant,length(x)/var_dim_constant)',1));
            new_asap_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * max(asap_power_sum_profile + [ones(1,prb.N_asap)*prb.station.pow_max zeros(1,var_dim_constant-prb.N_asap)]); 
            new_leave_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * max(asap_power_sum_profile)], v) + par.mu * (lse(z) - z' * prb.THETA' * v);
%     J = dot([new_flex_obj; 
%         new_asap_obj; 
%         new_leave_obj], v) + par.mu * (lse(z) - z' * prb.THETA' * v);
end
end