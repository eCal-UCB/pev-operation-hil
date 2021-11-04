import numpy as np
import cvxpy as cp
import datetime
from scipy import interpolate
import math
import timeit
from scipy.optimize import minimize
import sys

class Parameters:
    def __init__(self, TOU = np.ones((96,)) ,sens_num_poles = np.arange(2, 18, 3), monte_num_sims = 1, sim_starttime = 7, sim_endtime = 22,
                 sim_isFixedEventSequence = False, sim_isFixedSeed = False, sim_num_Events = 50, Ts = 0.25,
                 base_Tarriff_overstay = 1.0, station_num_poles = 8, eff = 0.92, lam_x = 10, lam_z_c = 10, lam_z_uc = 10,
                 lam_h_c = 10, lam_h_uc = 10, mu = 1e4, soft_v_eta = 1e-2, opt_eps = 1e-4,VIS_DETAIL = True):
        self.sens_analysis_num_poles = sens_num_poles
        # par.sens_analysis.num_poles = [2, 3]
        # monte carlo
        self.monte_num_sims = monte_num_sims
        # each one day simulation
        self.sim_starttime = sim_starttime
        self.sim_endtime = sim_endtime
        self.sim_isFixedEventSequence = sim_isFixedEventSequence
        self.sim_isFixedSeed = sim_isFixedSeed
        self.sim_num_events = sim_num_Events
        self.Ts = Ts
        # timestep, hour - - must decompose   1
        #  baseline  parameters
        self.base_tariff_overstay = base_Tarriff_overstay
        #TOU
        self.TOU = TOU

        # charging station config
        self.station_num_poles = station_num_poles
        # number  of  charging  poles
        self.eff = eff # power efficiency

        # dcm params
        self.dcm_choices = ['charging with flexibility', 'charging asap', 'leaving without charging']

        #pdfs
        self.pdf_visit = np.hstack((0.1 * np.ones((1, 7)), # 0 - 7
                         0.3 * np.ones((1, 5)), # 7 - 12
                         0.2 * np.ones((1, 2)), # 12 - 14
                         0.2 * np.ones((1, 2)), # 14 - 16
                         0.2 * np.ones((1, 6)), # 16 - 22
                         0.001 * np.ones((1, 2)))) # 22 - 24

        # regularization params
        self.lam_x = lam_x
        self.lam_z_c = lam_z_c
        self.lam_z_uc = lam_z_uc
        self.lam_h_c = lam_h_c
        # TODO: should be average overstay penalty in real data, should move to par
        self.lam_h_uc = lam_h_uc
        # TODO: should be average overstay penalty in real data, should move to par
        self.mu = mu
        self.soft_v_eta = soft_v_eta #softening equality constraint for v; to avoid numerical error
        self.opt_eps = opt_eps
        # debug_mode
        self.VIS_DETAIL = VIS_DETAIL

class Problem:
    def __init__(self, par ,**kwargs):
        self.par = par
        # if nargin == 0:
        #     # user input
        #     event = {}
        #     event["time"] = 14.25
        #     event["SOC_init"] = 0.3
        #     event["SOC_need"]= 0.5
        #     event["batt_cap"] = 80
        #     event["duration"] = 8
        #     event["overstay_duration"] = 1
        #     event["pow_max"] = 7.2
        #     event["pow_min"] = 0 
        # else:
        event = kwargs["event"]
        self.event = event
        
        self.user_time = round(event["time"] / par.Ts) * par.Ts
        self.user_SOC_init = event["SOC_init"]
        self.user_SOC_need = event["SOC_need"]
        self.user_batt_cap = event["batt_cap"]
        self.user_duration = round(event["duration"] / par.Ts) * par.Ts
        self.user_overstay_duration = round(event["overstay_duration"] / par.Ts) * par.Ts
        self.station_pow_max = event["pow_max"]
        self.station_pow_min = event["pow_min"]
        #dcm params
        # asc_flex = 2 + 0.2 * prb.user.duration;
        # asc_asap = 2.5;
        asc_flex = 2 + 0.401 * self.user_duration - 1.8531 * self.user_SOC_init #5.0583
        asc_asap = 1 + 0.865 * self.user_duration - 1.8531 * self.user_SOC_init #3.7088
        asc_leaving = 0
        energy_need = self.user_SOC_need * self.user_batt_cap
        # self.dcm_charging_flex.params = [-1 0 0 asc_flex].T
        # DCM parameters for choice 1 -- charging with flexibility
        # self.dcm_charging_asap.params = [0 - 1 0 asc_asap]';
        # DCM parameters for choice 2 -- charging as soon as possible
        self.dcm_charging_flex_params = np.array([[-0.1881 * energy_need], [0], [0], [asc_flex]])
        #% DCM parameters for choice 1 -- charging with flexibility
        self.dcm_charging_asap_params = np.array([[0], [- 0.1835 * energy_need], [0],[asc_asap]])
        #% DCM parameters for choice 2 -- charging as soon as possible
        self.dcm_leaving_params = np.array([[-0.01], [-0.01], [0.005], [asc_leaving]])
        #% DCM parameters for choice 3 -- leaving without charging
        self.THETA = np.vstack((self.dcm_charging_flex_params.T, self.dcm_charging_asap_params.T,
                     self.dcm_leaving_params.T))
        # problem specifications
        self.N_flex = int(self.user_duration / par.Ts) # charging duration that is not charged, hour
        
        ### IS THIS CORRECT? WHATS SOC NEED REPRESENTS? 
        self.N_asap = math.floor((self.user_SOC_need - self.user_SOC_init) *
                                 self.user_batt_cap / self.station_pow_max / par.eff / par.Ts)
        
        self.TOU = interpolate.interp1d(np.arange(0, 24 - 0.25 + 0.1, 0.25), par.TOU, kind = 'nearest')(np.arange(self.user_time,0.1 + self.user_time + self.user_duration - par.Ts,par.Ts)).T
class Optimization:

    def __init__(self, par, prb):
        self.Parameters = par
        self.Problem = prb
        self.opt_z = None
        self.opt_tariff_asap = None
        self.opt_tariff_flex = None
        self.opt_tariff_overstay = None

    def argmin_x(self, z, v):
        """
        Function to minimize charging cost. Flexible charging with variable power schedule
        Inputs: 

        Parameters: 
        z, array where [tariff_flex, tariff_asap, tariff_overstay, leave = 1 ]
        v, array with softmax results [sm_c, sm_uc, sm_y] (sm_y = leave)
        lam_x, regularization parameter for sum squares of the power var (u)
        lam_h_c, regularization parameter for g_flex
        lam_h_uc, regularization parameter for g_asap
        N_flex: timesteps arrival to departure 
        N_asap: timesteps required when charging at full capacity

        Parameters: 
        Decision Variables:
        SOC: state of charge (%)
        u: power (kW)

        Objective Function:
        Note: delta_k is not in the objective function!! 
        Check if it could make difference since power is kW cost is per kWh 

        """
    
        ### Read parameters 
        N_flex = self.Problem.N_flex
        N_asap = self.Problem.N_asap
        TOU = self.Problem.TOU
        station_pow_max = self.Problem.station_pow_max
        lam_x = self.Parameters.lam_x
        lam_h_c = self.Parameters.lam_h_c
        lam_h_uc = self.Parameters.lam_h_uc
        user_SOC_init  =  self.Problem.user_SOC_init
        user_SOC_need = self.Problem.user_SOC_need
        delta_k = self.Parameters.Ts
        eff = self.Parameters.eff
        user_bat_cap = self.Problem.user_batt_cap  

        print(N_flex, N_asap, TOU)

        
        # N_flex = kwargs["N_flex"]
        # N_asap = kwargs["N_asap"]
        # TOU = kwargs["TOU"]
        # station_pow_max =  kwargs["station_pow_max"]
        # lam_x = kwargs["lam_x"]
        # lam_h_c = kwargs["lam_h_c"]
        # lam_h_uc = kwargs["lam_h_uc"]
        # user_SOC_init  =   kwargs["user_SOC_init"]
        # user_SOC_need =  kwargs["user_SOC_need"]
        # delta_k =  kwargs["delta_k"]
        # eff =  kwargs["eff"]
        # user_bat_cap =  kwargs["user_bat_cap"]

        if sum(v) < 0 | (np.sum(v) < 1 - self.Parameters.soft_v_eta) | (np.sum(v) > 1 + self.Parameters.soft_v_eta):
            raise ValueError('[ ERROR] invalid $v$')
        

        ### Decision Variables
        SOC = cp.Variable(shape = (N_flex + 1))
        u = cp.Variable(shape = (N_flex))

        ### Define objective function
        # Flex Charging 
        #f_flex = cp.multiply(u , (TOU[:N_flex] - z[0]).reshape((N_flex))) # + cp.sum_squares(u) * lam_x
        f_flex = u.T @ TOU + cp.sum_squares(u) * lam_x
        g_flex = lam_h_c * 1 / z[2] 
        
        J_1 =  v[0] * (f_flex)
        
        # ASAP Charging
        f_asap = cp.sum(station_pow_max * (TOU[:N_asap] - z[1]))
        g_asap = lam_h_uc * 1 / z[2] 
        J_2 =  v[1] * (f_asap + g_asap)
        # Leave
        J_3 = cp.sum(TOU[:N_asap])

        J =   u.T @ TOU


        ## Constraints 

        constraints =  [SOC[0] == user_SOC_init]
        constraints += [SOC[N_flex] >=  user_SOC_need]
        constraints += [SOC <= 1]
        constraints += [u >= 0]
        constraints += [u <= station_pow_max]

        for i in range(0,N_flex ): 
            constraints += [SOC[i + 1] == SOC[i] + (eff * delta_k * u[i]) / user_bat_cap]


        ## Solve 
        obj = cp.Minimize(J)
        prob = cp.Problem(obj, constraints)
        prob.solve(solver='GUROBI')
        print("u:",np.round(u.value,2 ),"SOC:",np.round(SOC.value, 2))

        return  u.value, SOC.value

    def run_opt(self):
        start = timeit.timeit()

        itermax = 10000
        count = 0
        improve = np.inf
        zk = np.array([10,10,10,1]).reshape(4,1)
        # [z_c, z_uc, y, 1];
        # xk = np.ones((2 * self.Problem.N_flex + 1, 1)) # [soc0, ..., socN, u0, ..., uNm1]; - multiple dimensions 1 +  # of FLEX
        vk = np.array([0.45, 0.45, 0.1]).reshape(3,1)                     # [sm_c, sm_uc, sm_y]
        Jk = np.zeros((itermax, 1))

        power_flex, SOC_flex = self.argmin_x(zk, vk)

        opt = {}
        opt["z"] = zk
        opt["tariff_flex"] = zk[0]
        opt["tariff_asap"] = zk[1]
        opt["tariff_overstay"] = zk[2]
        # opt["x"] = xk
        # update demand charge
        opt["peak_pow"] = max(power_flex)
        opt["flex_SOCs"] = SOC_flex
        opt["flex_powers"] = power_flex
        opt["asap_powers"] = np.ones((self.Problem.N_asap, 1)) * self.Problem.station_pow_max
        opt["v"] = vk
        opt["prob_flex"] = vk[0]
        opt["prob_asap"] = vk[1]
        opt["prob_leave"] = vk[2]
        opt["J"] = Jk[0:count-1]
        opt["num_iter"] = count
        opt["prb"] = self.Problem
        opt["par"] = self.Parameters
        opt["time_start"] = self.Problem.user_time
        opt["time_end_flex"] = self.Problem.user_time + self.Problem.user_duration
        opt["time_end_asap"] = self.Problem.user_time + self.Problem.N_asap * self.Parameters.Ts
        end = timeit.timeit()
        return opt

def main(new_event = False, time = None, pow_min = None, pow_max = None, overstay_duration = None, duration= None,batt_cap =None, SOC_need = None, SOC_init = None):
    
    TOU = np.ones((96,)) 
    TOU[40:60] = np.ones(shape = (20,)) * 5
    par = Parameters(TOU = TOU)
    if new_event:
        event = {"time" : time, "pow_min" : pow_min, "pow_max" : pow_max, "overstay_duration" : overstay_duration, "duration" : duration, "batt_cap" : batt_cap, "SOC_need" : SOC_need, "SOC_init" : SOC_init}
        prb = Problem(par = par, event = event)
    else:
        prb = Problem(par=par)
    opt = Optimization(par, prb)

    z = [50, 50, 10, 1]
    v = [0.45,0.45,0.1]
    opt = opt.run_opt()

    return

if __name__ == "__main__":
    main(new_event = True, time = 4.25, pow_max = 10, pow_min = 0, overstay_duration= 1, duration= 12, batt_cap= 100, SOC_init = 0.3, SOC_need= 0.8)


