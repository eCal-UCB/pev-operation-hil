function xk = argmin_x(z,~,v)
par = get_glob_par();
prb = get_glob_prb();
N = prb.N_flex;

if (v < 0) | (sum(v) < 1-par.soft_v_eta) | (sum(v) > 1+par.soft_v_eta)
    error('[ ERROR] invalid $v$');
end

% cost function
J = @(x) dot([sum((x(N+2:end).*(prb.TOU(1:N) - z(1))).^2) + par.lambda.h_c * 1/z(3);
            sum((par.station.pow_max*(prb.TOU(1:N) - z(2))).^2) + par.lambda.h_uc * 1/z(3);
            sum((par.station.pow_max*(prb.TOU(1:N) - z(2))).^2)],v); 

% inequality constraint
A1L = [zeros(1,N) -1]; A1R = zeros(1,N);
A = [A1L A1R]; 
b = -prb.user.SOC_need;

% lower bound - power min
lb1 = zeros(N+1,1); lb1(end) = prb.user.SOC_need;
lb2 = par.station.pow_min.*ones(N,1);
lb = [lb1;lb2];

% upper bound - power max
ub1 = ones(N+1,1); 
ub2 = par.station.pow_max.*ones(N,1);
ub = [ub1;ub2];

% equality constraints - system dynamics
C1L = [1 zeros(1,N)]; C1R = zeros(1,N); % initial soc
C2L = [diag(-1.*ones(1,N)) zeros(N,1)] + [zeros(N,1)  diag(ones(1,N))];
C2R = -diag(par.eff/prb.user.batt_cap*ones(1,N));
d1 = prb.user.SOC_init; 
d2 = zeros(N,1);

Aeq = [C1L C1R; C2L C2R];
beq = [d1; d2];

% solve optimization
options = optimoptions('fmincon','Display','off');
xk = fmincon(J,prb.x0,A,b,Aeq,beq,lb,ub,[],options);
end