function vk = argmin_v(z,x,~)
par = get_glob_par();
prb = get_glob_prb();

% lse conjugate
lse_conj = @(v) dot(v,log(v));

% cost function
% J = @(v) dot([sum((x(prb.N_flex+2:end).*(prb.TOU(1:prb.N_flex) - z(1))).^2) + par.lambda.h_c * 1/z(3);
%             sum((par.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2))).^2) + par.lambda.h_uc * 1/z(3);
%             1/3*sum((par.station.pow_max*(prb.TOU(1:prb.N_asap) - 0)).^2)],v) ...
%          + par.mu * (lse_conj(v) - v' * prb.THETA * z);
J = @(v) dot([(sum((x(prb.N_flex+2:end).*(prb.TOU(1:prb.N_flex) - z(1)))+par.lambda.x.*x(prb.N_flex+2:end))+par.lambda.z_c*z(1)^2) + par.lambda.h_c * 1/z(3);
            sum((prb.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2)))+par.lambda.z_uc*z(2)^2) + par.lambda.h_uc * 1/z(3);
            sum(prb.station.pow_max*(prb.TOU(1:prb.N_asap) - 0))],v) ...
         + par.mu * (lse_conj(v) - v' * prb.THETA * z);
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
end