function xk = argmin_x(z,~,v)
global par
Ts = par.Ts;
N = par.N_flex;

% cost function
J_c = @(x) dot([sum((x(N+2:end).*(par.TOU(1:N*Ts) - z(1))).^2) + par.h_c.lambda * 1/z(3);
            sum((par.pow_max*(par.TOU(1:N*Ts) - z(2))).^2) + par.h_uc.lambda * 1/z(3);
            sum((par.pow_max*(par.TOU(1:N*Ts) - z(2))).^2)],v); %h2

% inequality constraint
A1L = diag(ones(1,N+1)); A1R = zeros(N+1,N); A1L(end,end) = -1;
A2L = zeros(N,N+1); A2R = diag(ones(1,N));
A3L = zeros(N,N+1); A3R = diag(-ones(1,N));
b1 = ones(N+1,1); b1(end) = -par.SOC_need;
b2 = par.pow_max .* ones(N,1);
b3 = -par.pow_min .* ones(N,1);

A = [A1L A1R; A2L A2R; A3L A3R];
b = [b1; b2; b3];

% equality constraints
C1L = zeros(1,N+1); C1R = zeros(1,N); C1L(1) = 1;
C2L = [diag(-1.*ones(1,N)) zeros(N,1)] + [zeros(N,1)  diag(ones(1,N))];
C2R = diag(par.Ts*par.eff/par.batt_cap*ones(1,N));
d1 = par.SOC_init; d2 = zeros(N,1);

C = [C1L C1R; C2L C2R];
d = [d1; d2];

% solve optimization
options = optimoptions('fmincon');
xk = fmincon(J_c,par.x0,A,b,C,d,[],[],[],options);
end