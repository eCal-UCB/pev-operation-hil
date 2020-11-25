#  Single charger optimization
import get_glob_prb
import get_glob_par
import numpy as np 
import cvxpy as cp 
import sys 

def run_opt():
	par = get_glob_par
	prb = get_glob_prb
	def J(z,x,v):
	    dm = np.dot([(sum((x[prb.N_flex+2:]*(prb.TOU(1:prb.N_flex) - z[1])) + par.lamb.x * x[prb.N_flex+2:]) + par.lamb.z_c*z[1]^2) + par.lamb.h_c * 1/z[3];

 	     	sum((prb.station.pow_max*(prb.TOU[1:prb.N_asap] - z[2])) + par.lamb.z_uc*z[2]^2) + par.lamb.h_uc * 1/z[3]; 

 	       	sum(prb.station.pow_max*(prb.TOU[1:prb.N_asap] - 0))], [v])
 	   return dm 
    itermax = 1e4
    count = 0
    improve = inf 
    zk = [0, 0, 0, 1]
    xk = np.ones(2*prb.N_flex+1,1+1)
    vk = [0.45, 0.45, 0.1]
    Jk = np.zeros(itermax, 1+1)

    while count < itermax < && improve >= 0 && abs(improve) >= par.opt.eps: 
    	count = count + 1

    	z = zk 
    	x = xk
    	v = vk 
    	Jk(count) = J(zk, xk, vk)

    	prb.z0 = zk;
    	prb.x0 = xk;
    	prb.v0 = vk;

    	set_glob_prb(prb);

    	zk = argmin_z([],xk,vk);
    	xk = argmin_x(zk,[],vk);
    	vk = argmin_v(zk,xk,[]);

    	improve = Jk(count)-J(zk,xk,vk)
    
    	opt.z = z;
		opt.tariff.flex = z[1];
		opt.tariff.asap = z[2];
		opt.tariff.overstay = z[3];
		opt.x = x;
		# update demand charge
		opt.peak_pow = max(x(prb.N_flex+2:));
		opt.flex.SOCs = x(1:prb.N_flex+1);
		opt.flex.powers = x(prb.N_flex+2:);
		opt.asap.powers = ones(prb.N_asap,1)*prb.station.pow_max;
		opt.prob.flex = v[1];
		opt.prob.asap = v[2];
		opt.prob.leave = v[3];
		opt.J = Jk[1:count];
		opt.num_iter = count;
		opt.prb = prb;
		opt.par = par;

		opt.time.start = prb.user.time;
		opt.time.end_flex = prb.user.time + prb.user.duration;
		opt.time.end_asap = prb.user.time + prb.N_asap*par.Ts;


		if par.VIR_DETAIL:
			print('')
    

