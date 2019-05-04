function zk = argmin_z(~,x,v)
par = get_glob_par();
prb = get_glob_prb();

if (v < 0) | (sum(v) < 1-par.soft_v_eta) | (sum(v) > 1+par.soft_v_eta)
    error('[ ERROR] invalid $v$');
end

% LSE function
lse = @(z) log(exp(dot(par.THETA(1,:),z))+exp(dot(par.THETA(2,:),z))+exp(dot(par.THETA(3,:),z)));

% cost function
J = @(z) dot([sum((x(prb.N_flex+2:end).*(prb.TOU(1:prb.N_flex) - z(1))).^2) + par.lambda.h_c * 1/z(3);
            sum((par.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2))).^2) + par.lambda.h_uc * 1/z(3);
            sum((par.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2))).^2)],v) ...
         + par.mu * (lse(z) - z' * par.THETA' * v);
     
% inequality constratins
A = [1 -1 0 0]; b = 0;
     
% lower and upper bound
lb = eps .* ones(4,1);
ub = 3 .* max(prb.TOU) .* ones(4,1);

% equality constraint
Aeq = [0 0 0 1]; beq = 1;

% solve optimization
options = optimoptions('fmincon','Display','off');
zk = fmincon(J,prb.z0,A,b,Aeq,beq,lb,ub,[],options);
end