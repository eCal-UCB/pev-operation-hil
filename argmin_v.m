function vk = argmin_v(z,x,~)
par = get_glob_par();

% lse conjugate
lse_conj = @(v) dot(v,log(v));

% cost function
J = @(v) dot([sum((x(par.N_flex+2:end).*(par.TOU(1:par.N_flex) - z(1))).^2) + par.lambda.h_c * 1/z(3);
            sum((par.station.pow_max*(par.TOU(1:par.N_asap) - z(2))).^2) + par.lambda.h_uc * 1/z(3);
            sum((par.station.pow_max*(par.TOU(1:par.N_asap) - z(2))).^2)],v) ...
         + par.mu * (lse_conj(v) - v' * par.THETA * z);

% inequality constraints
A = diag(-ones(1,3)); b = zeros(3,1);

% soft equality constraints
Aeq = [-ones(1,3);ones(1,3)]; beq = [-(1-par.soft_v_eta);1+par.soft_v_eta];
     
% solve optimization
options = optimoptions('fmincon','Display','off');
vk = fmincon(J,par.v0,[A;Aeq],[b;beq],[],[],[],[],[],options);
end