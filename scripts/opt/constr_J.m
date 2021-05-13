function J = constr_J(par,prb,z,x,v,station,k,existing_user_info,var_dim_constant)
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
        user = station(user_keys{i-1});
        overstay_cost = (user.time.leave - user.time.end) * user.z(3);
%         existing_flex_obj = existing_flex_obj + (sum(x(adj_constant+duration+2:adj_constant+2*duration+1,1).*(user.prb.TOU(TOU_idx:end) - user.price)) - overstay_cost);
        existing_flex_obj = existing_flex_obj + (sum(x(adj_constant+N_max+2:adj_constant+N_max+2+duration-1,1).*(user.prb.TOU(TOU_idx:end) - user.price)) - overstay_cost);
    end
    % existing asap user
    user_keys = station('ASAP_list');
    existing_asap_obj = 0;
    for i = 1:length(user_keys) % sum of users
        user = station(user_keys{i});
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
            if t <= opt.time.end % 04.21.2021 marked changed, correct from "k" to "t"
                if length(opt.powers) > 1
                    pow = interp1(opt.time.start:par.Ts:opt.time.end-par.Ts,opt.powers,t);
                else
                    pow = opt.powers;
                end
                try
                    asap_power_sum_profile(it) = asap_power_sum_profile(it) + pow;
                catch
                    a = 1;
                end
            end
        end
    end
    
    % part 1: case 1 - charging-FLEX
%     new_flex_obj = (sum((x(prb.N_flex+2:2*prb.N_flex+1,1).*(prb.TOU(1:prb.N_flex) - z(1)))...
%                         +par.lambda.x .* x(prb.N_flex+2:2*prb.N_flex+1,1))...
%                     +par.lambda.z_c*z(1)^2)...
%                   +par.lambda.h_c * 1/z(3); % with convergence regularization
    
    new_flex_obj = sum((x(N_max+2:N_max+2+prb.N_flex-1).*(prb.TOU(1:prb.N_flex) - z(1)))...
                +par.lambda.z_c*z(1)^2)...
              +par.lambda.h_c * 1/z(3); % with convergence regularization

    new_asap_obj = (sum((prb.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2)))...
                    +par.lambda.z_c*z(2)^2))...
                  +par.lambda.h_c * 1/z(3); % with convergence regularization

       
    % part 3: case 3 - leave
    new_leave_obj = sum(prb.station.pow_max*(prb.TOU(1:prb.N_asap) - 0));
    
    % part 4: demand charge
    all_power_profile = asap_power_sum_profile + sum(reshape(x,var_dim_constant,length(x)/var_dim_constant)',1);
    existing_power_profile = asap_power_sum_profile + sum(reshape(x(var_dim_constant+1:end),var_dim_constant,(length(x)/var_dim_constant)-1)',1);
    
    % part 5: power trajectory
    current_hour = floor(k) + 1; % instead of ceil(k)
    current_hour_energy = par.optSol.E_DA(current_hour) * par.optSol.to_MW * 1000;
    current_subhour_steps = (current_hour - k) / par.Ts + 1;
    
    % final objective function
    try
        J = dot([new_flex_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * max(all_power_profile);
            new_asap_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * max(existing_power_profile + [ones(1,prb.N_asap)*prb.station.pow_max zeros(1,var_dim_constant-prb.N_asap)]); 
            new_leave_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * max(asap_power_sum_profile); ...
            par.TOU_RT(k/par.Ts) * max(current_hour_energy - sum(all_power_profile(:,N_max+2:N_max+1+current_subhour_steps))*par.Ts, 0)], v);
    catch
        J = dot([new_flex_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * max(all_power_profile);
            new_asap_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * max(existing_power_profile + [ones(1,prb.N_asap)*prb.station.pow_max zeros(1,var_dim_constant-prb.N_asap)]); 
            new_leave_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * max(asap_power_sum_profile); ...
            par.TOU_RT(k/par.Ts) * max(current_hour_energy - sum(all_power_profile(:,N_max+2:end))*par.Ts, 0)], v);
    end
end

%J = @(z,x,v) dot([(sum((x(prb.N_flex+2:end).*(prb.TOU(1:prb.N_flex) - z(1)))+par.lambda.x.*x(prb.N_flex+2:end))+par.lambda.z_c*z(1)^2) + par.lambda.h_c * 1/z(3);
%            sum((prb.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2))) + par.lambda.z_uc*z(2)^2) + par.lambda.h_uc * 1/z(3);
%            sum(prb.station.pow_max*(prb.TOU(1:prb.N_asap) - 0))],v); % h_l