function xk = argmin_x_station(z,~,v,station,k,existing_user_info,var_dim_constant)
% flex_user_info(i, :) = [start_time, end_time, N_remain, TOU_idx, SOC_need, SOC_now];

par = get_glob_par();
prb = get_glob_prb();
N = prb.N_flex;

if (v < 0) | (sum(v) < 1-par.soft_v_eta) | (sum(v) > 1+par.soft_v_eta)
    error('[ ERROR] invalid $v$');
end

% cost function
% J = @(x) dot([sum((x(N+2:end).*(prb.TOU(1:N) - z(1))).^2) + par.lambda.h_c * 1/z(3);
%             sum((par.station.pow_max*(prb.TOU(1:N) - z(2))).^2) + par.lambda.h_uc * 1/z(3);
%             1/3*sum((par.station.pow_max*(prb.TOU(1:N) - 0)).^2)],v); 
% J = @(x) dot([(sum((x(prb.N_flex+2:end).*(prb.TOU(1:prb.N_flex) - z(1)))+par.lambda.x.*x(prb.N_flex+2:end))+par.lambda.z_c*z(1)^2) + par.lambda.h_c * 1/z(3);
%             sum((prb.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2)))+par.lambda.z_uc*z(2)^2) + par.lambda.h_uc * 1/z(3);
%             sum(prb.station.pow_max*(prb.TOU(1:prb.N_asap) - 0))],v); 
%         
J = @(xk) constr_J(xk);
num_flex_user = size(existing_user_info,1);

% inequality constraint
A_ineq_L = []; b_ineq = zeros(num_flex_user,1);
% lower bound - power min
lb = zeros(var_dim_constant*num_flex_user, 1);
% upper bound - power max
ub = zeros(var_dim_constant*num_flex_user, 1);
% equality constraints - system dynamics
A_eq = []; b_eq = [];
for idx = 1:size(existing_user_info,1)
    % iterate through new user + existing flex user
   N = existing_user_info(idx,3); % get user N_flex, aka the remian duration
   SOC_need = existing_user_info(idx,5);
   SOC_init = existing_user_info(idx,6);
   % ==== inequality constraint
   AL_ineq = [zeros(1,N) -1]; AR_ineq = zeros(1,var_dim_constant-N-1); % temporary L and R block
%    A_ineq(idx,:) = [AL_ineq AR_ineq];
   A_ineq_L = blkdiag(A_ineq_L, [AL_ineq AR_ineq]);
   b_ineq(idx) = -SOC_need; % SOC_need
   
   % lower bound - power min
   lb1 = zeros(N+1,1); lb1(end) = SOC_need;
   lb2 = prb.station.pow_min.*ones(N,1);
%    lb3 = zeros(2*(var_dim_constant-N),1);
   lb3 = zeros(var_dim_constant-2*N-1,1);
   lb((idx-1)*var_dim_constant+1:idx*var_dim_constant,1) = [lb1;lb2;lb3];
   
   % upper bound - power max
   ub1 = ones(N+1,1); 
   ub2 = prb.station.pow_max.*ones(N,1);
%    ub3 = zeros(2*(var_dim_constant-N),1);
   ub3 = zeros(var_dim_constant-2*N-1,1);
   ub((idx-1)*var_dim_constant+1:idx*var_dim_constant,1) = [ub1;ub2;ub3];
   
   % equality constraints - system dynamics
   C1L = [1 zeros(1,N)]; C1R = zeros(1,var_dim_constant-N-1); % initial soc
   C2L = [diag(-1.*ones(1,N)) zeros(N,1)] + [zeros(N,1)  diag(ones(1,N))];
   C2R = [-diag(par.eff*par.Ts/prb.user.batt_cap*ones(1,N)) zeros(N,var_dim_constant-2*N-1)];
   d1 = SOC_init; 
   d2 = zeros(N,1);

   A_eq = blkdiag(A_eq, [C1L C1R; C2L C2R]);
   b_eq = [b_eq; [d1; d2]];
end
% demand charge
C = var_dim_constant-1; % \tau=2 to \tau=T
M_dc_L = zeros(2*(C),size(A_ineq_L,2));
M_dc_R = zeros(2*(C),var_dim_constant);

for j = 1:C
    for i = 1:size(existing_user_info,1)
        % iterate through new user + existing flex user
        N = existing_user_info(i,3); % get user N_flex, aka the remian duration
        if j < N
            M_dc_L(2*(j-1)+1, (N+1)+j+(i-1)*var_dim_constant) = 1; % (N+1) skips the soc
        end
    end
    M_dc_R(2*(j-1)+1,j+1) = -1;
    M_dc_R(2*(j-1)+2,j:j+1) = [1 -1];
end
A_ineq_R = zeros(size(A_ineq_L,1), size(M_dc_R,2));
A_ineq = [A_ineq_L, A_ineq_R; M_dc_L, M_dc_R];

b2 = zeros(2*(C), 1);
% b2(1) = station('D_init'); b2(2) = -station('D_init');
num_asap_user = length(station('ASAP_list'));
for j = 1:C
    b2(2*(j-1)+1) = -num_asap_user * prb.station.pow_max;
end
b_ineq = [b_ineq; b2];

% fix dimension
A_eq = [A_eq, zeros(size(A_eq,1), size(A_ineq,2)-size(A_eq,2))];
lb = [lb; station('D_init'); zeros(size(A_ineq,2)-length(lb)-1,1)];
ub = [ub; station('D_init'); ones(size(A_ineq,2)-length(ub)-1,1)*station('pow_cap')];
% solve optimization
options = optimoptions('fmincon','Display','off');
xk = fmincon(J,ones(size(A_ineq,2),1),A_ineq,b_ineq,A_eq,b_eq,lb,ub,[],options);

function J = constr_J(x)
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
        existing_flex_obj = existing_flex_obj + (sum(x(adj_constant+duration+2:adj_constant+2*duration+1,1).*(user.prb.TOU(TOU_idx:end) - user.price)) + overstay_cost);
    end
    % existing asap user
    user_keys = station('ASAP_list');
    existing_asap_obj = 0;
    for i = 1:length(user_keys) % sum of users
        user = station(user_keys{1,i});
        overstay_cost = (user.time.leave - user.time.end) * user.z(3);
        TOU_idx = (k-user.time.start)/par.Ts+1;
        existing_asap_obj = existing_asap_obj + (sum(user.asap.powers*(user.prb.TOU(TOU_idx:end) - user.price)) + overstay_cost);
    end
    
    % ==== missing demand charge ====
    
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

    new_asap_obj = sum(prb.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2)))...
                  +par.lambda.h_c * 1/z(3); % without convergence regularization
%               +par.lambda.x .* x(prb.N_asap+2:2*prb.N_asap+1,1)
              
    % part 3: case 3 - leave
    new_leave_obj = sum(prb.station.pow_max*(prb.TOU(1:prb.N_asap) - 0));
    
    J = dot([new_flex_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * x(end); 
        new_asap_obj+existing_flex_obj+existing_asap_obj+station('cost_dc') * x(end); 
        new_leave_obj], v);
end
end