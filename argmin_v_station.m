function vk = argmin_v_station(z,x,~,station,k,existing_user_info,var_dim_constant)
par = get_glob_par();
prb = get_glob_prb();

% lse conjugate
lse_conj = @(v) dot(v,log(v));

% cost function
% J = @(v) dot([sum((x(prb.N_flex+2:end).*(prb.TOU(1:prb.N_flex) - z(1))).^2) + par.lambda.h_c * 1/z(3);
%             sum((par.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2))).^2) + par.lambda.h_uc * 1/z(3);
%             1/3*sum((par.station.pow_max*(prb.TOU(1:prb.N_asap) - 0)).^2)],v) ...
%          + par.mu * (lse_conj(v) - v' * prb.THETA * z);
% J = @(v) dot([(sum((x(prb.N_flex+2:end).*(prb.TOU(1:prb.N_flex) - z(1)))+par.lambda.x.*x(prb.N_flex+2:end))+par.lambda.z_c*z(1)^2) + par.lambda.h_c * 1/z(3);
%             sum((prb.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2)))+par.lambda.z_uc*z(2)^2) + par.lambda.h_uc * 1/z(3);
%             sum(prb.station.pow_max*(prb.TOU(1:prb.N_asap) - 0))],v) ...
%          + par.mu * (lse_conj(v) - v' * prb.THETA * z);

J = @(vk) constr_J(vk);
% inequality constraints
A = diag(-ones(1,3)); b = zeros(3,1);

% lower and upper bounds
lb = zeros(3,1);
ub = [1 1 0.3]';

% soft equality constraints
Aeq = [-ones(1,3);ones(1,3)]; beq = [-(1-par.soft_v_eta);1+par.soft_v_eta];
     
% solve optimization
options = optimoptions('fmincon','Display','off');
vk = fmincon(J,prb.v0,[A;Aeq],[b;beq],[],[],lb,ub,[],options);

function J = constr_J(v)
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
        existing_asap_obj = existing_asap_obj + (sum(user.asap.powers*(user.prb.TOU(TOU_idx:end) - user.price)) - overstay_cost);
    end
    
    % part 1: case 1 - charging-FLEX
%     new_flex_obj = (sum((x(prb.N_flex+2:2*prb.N_flex+1,1).*(prb.TOU(1:prb.N_flex) - z(1)))...
%                         +par.lambda.x .* x(prb.N_flex+2:2*prb.N_flex+1,1))...
%                     +par.lambda.z_c*z(1)^2)...
%                   +par.lambda.h_c * 1/z(3); % with convergence regularization

    new_flex_obj = (sum((x(prb.N_flex+2:2*prb.N_flex+1,1).*(prb.TOU(1:prb.N_flex) - z(1)))...
                        +par.lambda.x .* x(prb.N_flex+2:2*prb.N_flex+1,1)))...
                  +par.lambda.h_c * 1/z(3); % without convergence regularization 
              
    % part 2: case 2 - charging-ASAP
%     new_asap_obj = (sum((prb.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2)))...
%                         +par.lambda.x .* x(prb.N_asap+2:2*prb.N_asap+1,1))...
%                     +par.lambda.z_c*z(2)^2)...
%                   +par.lambda.h_c * 1/z(3); % with convergence regularization

    new_asap_obj = (sum((prb.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2)))...
                        +par.lambda.x .* x(prb.N_asap+2:2*prb.N_asap+1,1)))...
                  +par.lambda.h_c * 1/z(3); % without convergence regularization
              
    % part 3: case 3 - leave
    new_leave_obj = sum(prb.station.pow_max*(prb.TOU(1:prb.N_asap) - 0));
    
    J = dot([new_flex_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * x(end); 
        new_asap_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * x(end); 
        new_leave_obj], v) + par.mu * (lse_conj(v) - v' * prb.THETA * z);
end
end