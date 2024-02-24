import numpy as np
import cvxpy as cp
import datetime
from scipy import interpolate
import math
import timeit
from scipy.optimize import minimize
import sys
from datetime import datetime
import json 
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

class Parameters:
    ## Here modify the TOU to support overnight charging 
    def __init__(self, 
                z0 = np.array([1,1,1,1]).reshape(4,1) ,
                v0 = np.array([0.3333, 0.3333, 0.3333]).reshape(3,1) ,  
                Ts = 0.25, # interval size hour 
                base_tarriff_overstay = 1.0, 
                eff = 1, 
                soft_v_eta = 1e-4, 
                opt_eps = 0.0001, 
                TOU = np.ones((96,))):  
        ## TOU Tariff in cents / kwh 
        # TOU x power rate = cents  / hour 
        self.v0 = v0
        self.z0 = z0
        self.Ts = Ts
        self.base_tariff_overstay = base_tarriff_overstay
        self.TOU = TOU
        self.eff = eff # power efficiency
        self.dcm_choices = ['charging with flexibility', 'charging asap', 'leaving without charging']
        self.soft_v_eta = soft_v_eta #softening equality constraint for v; to avoid numerical error
        self.opt_eps = opt_eps
        
        assert len(self.TOU) == int(24 / self.Ts), "Mismatch between TOU cost array size and discretization steps"



class Problem:
    """
    time, int, user interval 
    duration, int, number of charging intervals 
    """
    def __init__(self, par ,**kwargs):
        self.par = par
        event = kwargs["event"]
        self.event = event
        
        self.user_time = event["time"]
        self.e_need = event["e_need"]

        self.user_duration = event["duration"]
        # self.user_overstay_duration = round(event["overstay_duration"] / par.Ts) * par.Ts
        self.station_pow_max = event["station_pow_max"]
        # self.station_pow_min = event["pow_min"]
        self.user_power_rate = event['user_power_rate']

        self.power_rate = min(self.user_power_rate,self.station_pow_max)

        self.dcm_charging_flex_params = np.array([[ - self.power_rate * 0.0184 / 2], [ self.power_rate * 0.0184 / 2], [0], [0]])
        #% DCM parameters for choice 1 -- charging with flexibility
        self.dcm_charging_asap_params = np.array([[self.power_rate * 0.0184 / 2], [- self.power_rate * 0.0184 / 2], [0],[0.341 ]])
        #% DCM parameters for choice 2 -- charging as soon as possible
        self.dcm_leaving_params = np.array([[self.power_rate * 0.005 / 2], [self.power_rate * 0.005 / 2], [0], [-1 ]])


        
        #% DCM parameters for choice 3 -- leaving without charging
        self.THETA = np.vstack((self.dcm_charging_flex_params.T, self.dcm_charging_asap_params.T,
                     self.dcm_leaving_params.T))
        # problem specifications
        self.N_flex = self.user_duration # charging duration that is not charged, hour
        
        ### IS THIS CORRECT? WHATS SOC NEED REPRESENTS? 
        # self.N_asap = math.floor((self.user_SOC_need - self.user_SOC_init) *
        #                          self.user_batt_cap / self.station_pow_max / par.eff / par.Ts)

        ## HERE 12 IS SELF CODED 
        self.N_asap = math.ceil((self.e_need / self.power_rate / par.eff * int(1 / par.Ts)))
        self.N_asap_remainder = (self.e_need / self.power_rate / par.eff * int(1 / par.Ts)) % 1

#         print(par.TOU) 
        if len(par.TOU) < self.user_time + self.user_duration: # if there is overnight chaarging 
            

            par.TOU = np.concatenate([par.TOU,par.TOU]) 

        self.TOU = par.TOU[self.user_time:(self.user_time + self.user_duration)]
#         print(self.TOU) 
        # self.TOU = interpolate.interp1d(np.arange(0, 24 - 0.25 + 0.1, 0.25), par.TOU, kind = 'nearest')(np.arange(self.user_time, 0.1 + self.user_time + self.user_duration - par.Ts, par.Ts)).T
        
        assert self.N_asap <= self.N_flex, print("Not enought time (n_asap,n_flex)",self.N_asap,self.N_flex)

class Optimization:

    def __init__(self, par, prb):
        self.Parameters = par
        self.Problem = prb
        self.opt_z = None
        self.opt_tariff_asap = None
        self.opt_tariff_flex = None
        self.opt_tariff_overstay = None

    def argmin_v(self, u, z):

        """

        Parameters 
        Decision Variables: 
        v: price [ sm(theta_flex, z), sm(theta_asap, z), sm(theta_leave, z) ], (3,1)

        """
  
        ### Read parameters 



        N_flex = self.Problem.N_flex
        N_asap = self.Problem.N_asap
        N_asap_remainder  = self.Problem.N_asap_remainder 
        TOU = self.Problem.TOU
        station_pow_max = self.Problem.power_rate

        # mu = self.Parameters.mu
        THETA = self.Problem.THETA 
        soft_v_eta = self.Parameters.soft_v_eta
        delta_t = self.Parameters.Ts




        ### Decision Variables
        v = cp.Variable(shape = (3), pos = True)

        ### Define objective function
        # Flex Charging 
        # reg_flex =  cp.norm(u,2) * lam_x + cp.norm(z[0],2) * lam_z_c 
        
        f_flex = u.T @ (TOU - z[0]) * delta_t
            

        J_1 =  v[0] * (f_flex)
        
        # ASAP Charging
        # reg_asap =  cp.norm(z[1],2) * lam_z_uc 

        # We
        if N_asap_remainder > 0:
            if N_asap <= 1:
                f_asap = ((station_pow_max * N_asap_remainder) * (TOU[N_asap - 1] - z[1])) * delta_t 
            # print(self.Problem.N_asap)
            # print(N_asap_remainder, N_asap_remainder.shape)
            # print(TOU[:N_asap - 1], TOU[:N_asap - 1].shape)
            else:
                f_asap = (cp.sum(station_pow_max * (TOU[:N_asap - 1] - z[1])) + (station_pow_max * N_asap_remainder) * (TOU[N_asap - 1] - z[1]) )* delta_t 

        else: 
            f_asap = cp.sum(station_pow_max * (TOU[:N_asap ] - z[1])) * delta_t 


        J_2 =  v[1] * (f_asap )
        
        # Leave
        # J_3 = v[2] * cp.sum(TOU[:N_asap] * station_pow_max) * delta_t
        J_3 = 0
        J =    J_1 + J_2 + J_3 

        ### Log sum function conjugate: negative entropy 
        # lse_conj = - cp.sum(cp.entr(v))
        # func = v.T @ (THETA @ z)
        # # J_4 = mu * (lse_conj - func) 
        # constraints += [ v <= np.array((1,1,1))] # What is this? 
        # constraints += [ cp.sum(v) >= 1 - soft_v_eta ]

        constraints = [v >= 0 ]
        constraints += [cp.sum(v) == 1 ]
        # constraints += [v[2] <= 0.50 ]

        constraints += [ cp.log_sum_exp(THETA @ z) - cp.sum(cp.entr(v)) - v.T @ (THETA @ z) <= soft_v_eta ]
        
        ## Solve 
        obj = cp.Minimize(J)
        prob = cp.Problem(obj, constraints)
        prob.solve()  

        # try:
        #     # print(  "v",v.value)
        #     # print(  "status",prob.status)
        #     temp = v.value
        # except:
        #     print(  "status",prob.status)
        return np.round(v.value,4)

    def argmin_z(self, u, v):
        """
        Function to determine prices 

        Decision Variables: 
        z: price [tariff_flex, tariff_asap, tariff_overstay, leave = 1 ]

        Parameters: 
        u, array, power for flex charging 
        v, array with softmax results [sm_c, sm_uc, sm_y] (sm_y = leave)
        lam_x, regularization parameter for sum squares of the power var (u)
        lam_z_c, regularization parameter for sum squares of the price flex (u)
        lam_z_uc, regularization parameter for sum squares of the price asap (u)
        lam_h_c, regularization parameter for g_flex
        lam_h_uc, regularization parameter for g_asap
        N_flex: timesteps arrival to departure 
        N_asap: timesteps required when charging at full capacity

        """
        if sum(v) < 0 | (np.sum(v) < 1 - self.Parameters.soft_v_eta) | (np.sum(v) > 1 + self.Parameters.soft_v_eta):
            raise ValueError('[ ERROR] invalid $v$')
        
        ### Read parameters 
        N_flex = self.Problem.N_flex
        N_asap = self.Problem.N_asap
        TOU = self.Problem.TOU
        station_pow_max = self.Problem.power_rate
#         vehicle_power_rate = self.Problem.power_rate
        delta_t = self.Parameters.Ts
        soft_v_eta = self.Parameters.soft_v_eta
        THETA = self.Problem.THETA 

        ### Decision Variables
        
        z = cp.Variable(shape = (4), pos = True)

        f_flex = u.T @ (TOU - z[0]) * delta_t
        # g_flex = lam_h_c * cp.inv_pos(z[2])

        J_1 =  v[0] * (f_flex)
        
        # ASAP Charging
        N_asap_remainder  = self.Problem.N_asap_remainder 
        # print(N_asap_remainder)
        
        if N_asap_remainder > 0:
            if N_asap <= 1:
                f_asap = ((station_pow_max * N_asap_remainder) * (TOU[N_asap - 1] - z[1])) * delta_t 
            # print(self.Problem.N_asap)
            # print(N_asap_remainder, N_asap_remainder.shape)
            # print(TOU[:N_asap - 1], TOU[:N_asap - 1].shape)
            else:
                f_asap = (cp.sum(station_pow_max * (TOU[:N_asap - 1] - z[1])) + (station_pow_max * N_asap_remainder) * (TOU[N_asap - 1] - z[1]) )* delta_t 

        else: 
            f_asap = cp.sum(station_pow_max * (TOU[:N_asap ] - z[1])) * delta_t 

        # g_asap =  lam_h_c* cp.inv_pos(z[2])
        
        J_2 =  v[1] * (f_asap)
        # Leave
        J_3 = 0

        J =    J_1 + J_2 + J_3 


        ### Log sum function 
        # lse = cp.log_sum_exp(THETA @ z)
        # func = z.T @ (THETA.T @ v)
        # J_4 = mu * (lse - func) 
        constraints = [z[3] == 1]
        constraints += [ cp.log_sum_exp(THETA @ z) - cp.sum(cp.entr(v)) - v.T @ (THETA @ z) <= soft_v_eta ]
        
        ## Solve 
        obj = cp.Minimize(J)
        prob = cp.Problem(obj, constraints)

        prob.solve()  

        # try:
        #     # print("z",np.round(z.value,5))
        #     temp = np.round(z.value,5)
        # except:
        #     print(  "z status",prob.status)
        
        
        return z.value


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

        Outputs
        u: power 
        SOC: SOC level 
        """
    
        ### Read parameters 
        N_flex = self.Problem.N_flex
        N_asap = self.Problem.N_asap
        TOU = self.Problem.TOU
        station_pow_max = self.Problem.power_rate
        vehicle_power_rate = self.Problem.power_rate
        e_need = self.Problem.e_need
        eff = 1
        # user_bat_cap = self.Problem.user_batt_cap  
        delta_t = self.Parameters.Ts 
        soft_v_eta = self.Parameters.soft_v_eta

        if sum(v) < 0 | (np.sum(v) < 1 - self.Parameters.soft_v_eta) | (np.sum(v) > 1 + self.Parameters.soft_v_eta):
            raise ValueError('[ ERROR] invalid $v$')
        

        ### Decision Variables
        e_delivered = cp.Variable(shape = (N_flex + 1))
        u = cp.Variable(shape = (N_flex))
        f_flex = u.T @ (TOU - z[0]) * delta_t
        # g_flex = lam_h_c * cp.inv_pos(z[2])
        
        J_1 =  v[0] * (f_flex)
        
        # ASAP Charging
        # reg_asap =  cp.norm(z[1],2) * lam_z_uc 
        # f_asap = cp.sum(station_pow_max * (TOU[:N_asap] - z[1])) * delta_t
        N_asap_remainder  = self.Problem.N_asap_remainder 
   
        if N_asap_remainder > 0:
            if N_asap <= 1:
                f_asap = ((station_pow_max * N_asap_remainder) * (TOU[N_asap - 1] - z[1])) * delta_t 
            # print(self.Problem.N_asap)
            # print(N_asap_remainder, N_asap_remainder.shape)
            # print(TOU[:N_asap - 1], TOU[:N_asap - 1].shape)
            else:
                f_asap = (cp.sum(station_pow_max * (TOU[:N_asap - 1] - z[1])) + (station_pow_max * N_asap_remainder) * (TOU[N_asap - 1] - z[1]) )* delta_t 

        else: 
            f_asap = cp.sum(station_pow_max * (TOU[:N_asap ] - z[1])) * delta_t 

        # g_asap = lam_h_uc * cp.inv_pos(z[2])
        J_2 =  v[1] * (f_asap )
        # Leave
        # J_3 = v[2] * cp.sum(TOU[:N_asap] * station_pow_max * delta_t) 
        J_3 = 0

        J =    J_1 + J_2 + J_3


        ## Constraints 

        constraints =  [e_delivered[0] == 0]
        constraints += [e_delivered[N_flex] >=  e_need]
        constraints += [e_delivered <= e_need]
        constraints += [u >= 0]
        constraints += [u <= station_pow_max ]

        # System dynamics
        for i in range(0,N_flex ): 
            constraints += [e_delivered[i + 1] == e_delivered[i] + (eff * delta_t * u[i])]


        ## Solve 
        obj = cp.Minimize(J)
        prob = cp.Problem(obj, constraints)
        prob.solve()
        
        # try:
        #     print("u:",np.round(u.value,2 ))
        # except:
        #     print(  "status",prob.status)
        return  u.value, e_delivered.value

    def run_opt(self):
        start = timeit.timeit()
        
        def J_func(z, u, v):
            ### Read parameters 
            N_asap = self.Problem.N_asap
            TOU = self.Problem.TOU
            station_pow_max = self.Problem.power_rate
            delta_t = self.Parameters.Ts 
            soft_v_eta = self.Parameters.soft_v_eta

            # reg_flex =  np.linalg.norm(u,2) * lam_x + z[0]**2 * lam_z_c
 
            f_flex = u.T @ (TOU - z[0]) * delta_t
            # g_flex = lam_h_c * 1 / z[2] 
            
            J_1 =  v[0] * (f_flex)
            
            # ASAP Charging
            # reg_asap =  z[1]**2 * lam_z_uc 
            # f_asap = np.sum(station_pow_max * (TOU[:N_asap] - z[1])) * delta_t
            N_asap_remainder  = self.Problem.N_asap_remainder 
            
            if N_asap_remainder > 0:
                f_asap = (np.sum(station_pow_max * (TOU[:N_asap - 1] - z[1])) + (station_pow_max * N_asap_remainder) * (TOU[N_asap - 1] - z[1]) )* delta_t 
            else: 
                f_asap = np.sum(station_pow_max * (TOU[:N_asap ] - z[1])) * delta_t 
    

            # g_asap = lam_h_uc * 1 / z[2] 
            J_2 =  v[1] * (f_asap )
            
            # Leave
            # Include the p_max 
            # J_3 = v[2] * np.sum(TOU[:N_asap]) * station_pow_max * delta_t
            J_3 = 0

            return  np.array([J_1 , J_2 , J_3])
        
        def charging_revenue(z, u):
            
            N_asap = self.Problem.N_asap
            TOU = self.Problem.TOU
            station_pow_max = self.Problem.power_rate

            delta_t = self.Parameters.Ts 

            f_flex = u.T @ (z[0]- TOU) * delta_t
            ## u : kW , z: cents / kWh, TOU : cents / kWh , delta_t : 1 \ h
            # f_asap = np.sum(station_pow_max * (z[1] - TOU[:N_asap])) * delta_t 

            N_asap_remainder  = self.Problem.N_asap_remainder 
            
            if N_asap_remainder > 0:
                f_asap = (np.sum(station_pow_max * (TOU[:N_asap - 1] - z[1])) + (station_pow_max * N_asap_remainder) * (TOU[N_asap - 1] - z[1]) )* delta_t 
                
            else: 
                f_asap = np.sum(station_pow_max * (TOU[:N_asap ] - z[1])) * delta_t 
                


            return f_flex, f_asap

        itermax = 1000
        count = 0
        improve = np.inf
        
        # [z_c, z_uc, y, 1];
        # xk = np.ones((2 * self.Problem.N_flex + 1, 1)) # [soc0, ..., socN, u0, ..., uNm1]; - multiple dimensions 1 +  # of FLEX
        
        zk = self.Parameters.z0
        uk_flex = np.zeros((self.Problem.N_flex))
        vk = self.Parameters.v0                  # [sm_c, sm_uc, sm_y]

        ###     THIS VALUES ARE STORED FOR DEBUGGING     ## 
        
        Jk = np.zeros((itermax))
        rev_flex = np.zeros((itermax))
        rev_asap = np.zeros((itermax))
        z_iter = np.zeros((4,itermax))
        v_iter = np.zeros((3,itermax))
        J_sub = np.zeros((3,itermax))


        # print(J_func(zk, uk_flex, vk))

        while (count < itermax) & (improve >= 0) & (abs(improve) >= 0.00001):

            Jk[count]  = J_func(zk, uk_flex, vk).sum()
            J_sub[:,count] = J_func(zk, uk_flex, vk).reshape(3,)    
            rev_flex[count], rev_asap[count] = charging_revenue(zk, uk_flex)
            z_iter[:,count] = zk.reshape((4,))
            v_iter[:,count] = vk.reshape((3,))

            uk_flex, e_deliveredk_flex = self.argmin_x(zk, vk)
            
            vk = self.argmin_v(uk_flex, zk)
            zk = self.argmin_z(uk_flex, vk)

            # compute residual
            # print(Jk[count])

            improve = Jk[count] - J_func(zk, uk_flex, vk).sum()
            # print(J_func(zk, uk_flex, vk))
            count += 1

        print(count, improve)
        opt = {}
        opt['e_need'] = self.Problem.e_need
        opt["z"] = zk
        opt["z_hr"] = zk * self.Problem.power_rate
        opt["tariff_flex"] = zk[0]
        opt["tariff_asap"] = zk[1]
        opt["tariff_overstay"] = zk[2]
        # opt["x"] = xk
        # update demand charge
        opt["peak_pow"] = max(uk_flex)
        opt["flex_e_delivered"] = e_deliveredk_flex
        opt["flex_powers"] = uk_flex
        asap_powers = np.ones((self.Problem.N_asap, 1)) * self.Problem.station_pow_max
        if self.Problem.N_asap_remainder != 0:
            asap_powers[self.Problem.N_asap - 1] =  self.Problem.station_pow_max * self.Problem.N_asap_remainder
        opt["asap_powers"] = asap_powers
        opt["v"] = vk
        opt["prob_flex"] = vk[0]
        opt["prob_asap"] = vk[1]
        opt["prob_leave"] = vk[2]

        opt["J"] = Jk[:count]
        opt["J_sub"] = J_sub[:,:count]
        opt["z_iter"] = z_iter[:,:count]
        opt["v_iter"] = v_iter[:,:count]
        opt["rev_flex"] =rev_flex[:count]
        opt["rev_asap"] = rev_asap[:count]

        opt["num_iter"] = count
        opt["prb"] = self.Problem
        opt["par"] = self.Parameters
        opt["time_start"] = self.Problem.user_time
        opt["time_end_flex"] = self.Problem.user_time + self.Problem.user_duration
        opt["time_end_asap"] = self.Problem.user_time + self.Problem.N_asap 
        end = timeit.timeit()
        return opt

def save_look_up_tables(output_folder_name, reg_price, sch_price, sch_power):
    """ Outputs 6 files """
    
    # check if outputs folder exists 
    MYDIR = ("Outputs")
    CHECK_FOLDER = os.path.isdir(MYDIR)

    # If folder doesn't exist, then create it.
    if not CHECK_FOLDER:
        os.makedirs(MYDIR)
        print("created folder : ", MYDIR)
    
    MYDIR = ("Outputs/Recent")
    CHECK_FOLDER = os.path.isdir(MYDIR)

    if not CHECK_FOLDER:
        os.makedirs(MYDIR)
        print("created folder : ", MYDIR)

        reg_price.to_excel( os.path.join("Outputs/Recent","Regular_Prices.xlsx") )
        sch_price.to_excel( os.path.join("Outputs/Recent","Scheduled_Prices.xlsx") )

    # check if outputs folder exists 
    MYDIR = (os.path.join("Outputs",output_folder_name))
    CHECK_FOLDER = os.path.isdir(MYDIR)

    # If folder doesn't exist, then create it.
    if not CHECK_FOLDER:
        os.makedirs(MYDIR)
        print("created folder : ", MYDIR)
        
        reg_price.to_excel( os.path.join(MYDIR,"Regular_Prices.xlsx") )
        sch_price.to_excel( os.path.join(MYDIR,"Scheduled_Prices.xlsx") )
        sch_power.to_excel( os.path.join(MYDIR,"Scheduled_Powers.xlsx") )

    elif CHECK_FOLDER:
        print("already exists: ", MYDIR)
        reg_price.to_excel( os.path.join(MYDIR,"Regular_Prices_1.xlsx") )
        sch_price.to_excel( os.path.join(MYDIR,"Scheduled_Prices_1.xlsx") )
        sch_power.to_excel( os.path.join(MYDIR,"Scheduled_Powers_1.xlsx") )



def main():

    ## This is the function that will run in every 30 minutes 
    z0 = np.array([25,35,1,1]).reshape(4,1) 

    ################## Define the TOU Cost ##################
    # off-peak 0.175  cents / kwh 
    TOU_tariff = np.ones((96,)) * 17.5
    ## 4 pm - 9 pm peak 0.367 cents / kwh 
    TOU_tariff[64:84] = 36.7
    ## 9 am - 2 pm super off-peak 0.49 $ / kWh  to cents / kwh
    TOU_tariff[36:56] = 14.9

    # Define the current time and the arrival interval 
    current_time = datetime.now()
    current_time_str = current_time.strftime("%H:%M")
    
    arrival_hour = int(current_time_str.split(':')[0])
    arrival_minute = int(np.floor(int(current_time_str.split(':')[1]) / 30) * 30)

    intervals_in_hour = 4 
    minutes_in_interval = 15 
    ### If the EV arrives sometime between arrival interval start and until next interval we will grab this information 
    arrival_interval = int(arrival_hour  * intervals_in_hour  + arrival_minute / minutes_in_interval )

    duration_list = np.arange(1,41) # for each optimization run for 10 hrs ahead with 15 min intervals 
    energy_list = np.arange(1,60,5) 
    # for each optimization run for 10 hrs ahead with 15 min intervals 


    ##### Create Output DataFrames #####
    output_folder_name = current_time.strftime("%m_%d_%y_") + str(arrival_hour) + "_" +  str(arrival_minute).ljust(2, '0')


    reg_price = pd.DataFrame(index = duration_list, columns =  energy_list, dtype=float)
    sch_price = pd.DataFrame(index = duration_list, columns =  energy_list, dtype=float)
    sch_power = pd.DataFrame(index = duration_list, columns =  energy_list, dtype=str)
    # for each optimization run for 10 hrs ahead with 15 min intervals 

    power_rate = 6.6

    for duration_interval in duration_list: 
        for e_needed in energy_list:
            if duration_interval < int(np.ceil((e_needed / power_rate) * intervals_in_hour)): # Timeste dependent value 
                print("Not enough time")
                continue
            else: 
                event = {
                        "time" : arrival_interval, 
                        "pow_min" : 0, 
                        "pow_max" : 6.6, 
                        "overstay_duration":1,
                        "duration" : duration_interval, 
                        "power_rate":power_rate, 
                        "e_need":e_needed 
                        }
                par = Parameters(  TOU = TOU_tariff, 
                                        z0 = z0,
                                        eff = 1.0, 
                                        Ts = 0.25) 
                
                prb = Problem(par = par, event = event)
                
                opt = Optimization(par, prb)
                
                res = opt .run_opt()

                reg_price.loc[duration_interval, e_needed] = res['tariff_asap'] * power_rate # tariff_asap in cents / kwh we need to convert cents / hour 
                sch_price.loc[duration_interval, e_needed] = res['tariff_flex'] * power_rate # tariff_asap in cents / kwh we need to convert cents / hour 
                sch_power.loc[duration_interval, e_needed] = json.dumps(list(res['flex_powers']))
                
    save_look_up_tables(output_folder_name, reg_price, sch_price, sch_power)

    return

if __name__ == "__main__":
    main()


