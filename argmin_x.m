function xk = argmin_x(z,~,v)
par = get_glob_par();
N = par.N_flex;

% cost function
J = @(x) dot([sum((x(N+2:end).*(par.TOU(1:N) - z(1))).^2) + par.lambda.h_c * 1/z(3)^2;
            sum((par.station.pow_max*(par.TOU(1:N) - z(2))).^2) + par.lambda.h_uc * 1/z(3)^2;
            sum((par.station.pow_max*(par.TOU(1:N) - z(2))).^2)],v); %h2

% inequality constraint
% A1L = diag(ones(1,N+1)); A1R = zeros(N+1,N); A1L(end,end) = -1;
% A2L = zeros(N,N+1); A2R = diag(ones(1,N));
% A3L = zeros(N,N+1); A3R = diag(-ones(1,N));
% b1 = ones(N+1,1); b1(end) = -par.user.SOC_need;
% b2 = par.station.pow_max .* ones(N,1);
% b3 = -par.station.pow_min .* ones(N,1);
% A = [A1L A1R; A2L A2R; A3L A3R];
% b = [b1; b2; b3];

A1L = [zeros(1,N) -1]; A1R = zeros(1,N);
A = [A1L A1R]; 
b = -par.user.SOC_need;

% lower bound - power min
lb1 = zeros(N+1,1); lb1(end) = par.user.SOC_need;
lb2 = par.station.pow_min.*ones(N,1);
lb = [lb1;lb2];

% upper bound - power max
ub1 = ones(N+1,1); 
ub2 = par.station.pow_max.*ones(N,1);
ub = [ub1;ub2];

% equality constraints - system dynamics
C1L = [1 zeros(1,N)]; C1R = zeros(1,N); % initial soc
C2L = [diag(-1.*ones(1,N)) zeros(N,1)] + [zeros(N,1)  diag(ones(1,N))];
C2R = -diag(par.eff/par.user.batt_cap*ones(1,N));
d1 = par.user.SOC_init; 
d2 = zeros(N,1);

Aeq = [C1L C1R; C2L C2R];
beq = [d1; d2];

% solve optimization
options = optimoptions('fmincon','Display','off');
xk = fmincon(J,par.x0,A,b,Aeq,beq,lb,ub,[],options);
end