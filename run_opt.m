function opt = run_opt()
tic;
par = get_glob_par();
prb = get_glob_prb();

%% Functions
% fprintf('TEST: N_asap: %d\n',prb.N_asap);
% fprintf('TEST: length(TOU): %d\n',length(prb.TOU));
J = @(z,x,v) dot([(sum((x(prb.N_flex+2:end).*(prb.TOU(1:prb.N_flex) - z(1)))+par.lambda.x.*x(prb.N_flex+2:end))+par.lambda.z_c*z(1)^2) + par.lambda.h_c * 1/z(3);
            sum((prb.station.pow_max*(prb.TOU(1:prb.N_asap) - z(2))) + par.lambda.z_uc*z(2)^2) + par.lambda.h_uc * 1/z(3);
            sum(prb.station.pow_max*(prb.TOU(1:prb.N_asap) - 0))],v); % h_l
%% Run algorithm -- block coordinate descent
itermax = 1e4;
count = 0; improve = inf;
zk = ones(4,1);                         % [z_c, z_uc, y, 1];
xk = ones(2*prb.N_flex+1,1);            % [soc0, ..., socN, u0, ..., uNm1]; - multiple dimensions 1+#of FLEX
vk = 1/3*ones(3,1);                     % [sm_c, sm_uc, sm_y];
Jk = zeros(itermax,1);
while count < itermax && improve >= 0 && abs(improve) >= par.opt.eps
    count = count + 1;
    Jk(count) = J(zk,xk,vk);    
    
    % update init variables
    prb.z0 = zk; prb.x0 = xk; prb.v0 = vk; set_glob_prb(prb);
    
    % update control variables
    zk = argmin_z([],xk,vk);
    xk = argmin_x(zk,[],vk);
    vk = argmin_v(zk,xk,[]);
    
    % compute residual
    improve = Jk(count)-J(zk,xk,vk);
    
%     if mod(count,1) == 0
%         fprintf('[ OPT] iter: %d, improve: %.3f\n',count,improve);
%     end
end

opt.z = zk;
opt.tariff.flex = zk(1);
opt.tariff.asap = zk(2);
opt.tariff.overstay = zk(3);
opt.x = xk;
% update demand charge
opt.peak_pow = max(xk(prb.N_flex+2:end));
opt.flex.SOCs = xk(1:prb.N_flex+1);
opt.flex.powers = xk(prb.N_flex+2:end);
opt.asap.powers = prb.station.pow_max;
opt.v = vk;
opt.prob.flex = vk(1);
opt.prob.asap = vk(2);
opt.prob.leave = vk(3);
opt.J = Jk(1:count);
opt.num_iter = count;
opt.prb = prb;
opt.par = par;

opt.time.start = prb.user.time;
opt.time.end_flex = prb.user.time + prb.user.duration;
opt.time.end_asap = prb.user.time + prb.N_asap*par.Ts;

fprintf('[%s OPT] DONE (%.2f sec) sum(vk) = %.2f, iterations = %d\n',datetime('now'),toc,sum(vk),count);
end