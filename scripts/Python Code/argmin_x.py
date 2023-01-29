import cvxpy as cp
# from opt import main
import numpy as np 
import matplotlib.pyplot as plt
import os 

def argmin_x(z, v, kwargs):

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
    # N_flex = self.Problem.N_flex
    # N_asap = self.Problem.N_asap
    # TOU = self.Problem.TOU
    # station_pow_max = self.Problem.station_pow_max
    # lam_x = self.Parameters.lam_x
    # lam_h_c = self.Parameters.lam_h_c
    # lam_h_uc = self.Parameters.lam_h_uc
    # user_SOC_init  =  self.Problem.user_SOC_init
    # user_SOC_need = self.Problem.user_SOC_need
    # delta_k = self.Parameters.Ts
    # eff = self.Parameters.eff
    # user_bat_cap = self.Problem.user_batt_cap  

    
    N_flex = kwargs["N_flex"]
    N_asap = kwargs["N_asap"]
    TOU = kwargs["TOU"]
    station_pow_max =  kwargs["station_pow_max"]
    lam_x = kwargs["lam_x"]
    lam_h_c = kwargs["lam_h_c"]
    lam_h_uc = kwargs["lam_h_uc"]
    user_SOC_init  =   kwargs["user_SOC_init"]
    user_SOC_need =  kwargs["user_SOC_need"]
    delta_k =  kwargs["delta_k"]
    eff =  kwargs["eff"]
    user_bat_cap =  kwargs["user_bat_cap"]

    # if sum(v) < 0 | (np.sum(v) < 1 - self.Parameters.soft_v_eta) | (np.sum(v) > 1 + self.Parameters.soft_v_eta):
    #     raise ValueError('[ ERROR] invalid $v$')
    

    ### Decision Variables
    SOC = cp.Variable(shape = (N_flex + 1))
    u = cp.Variable(shape = (N_flex))

    ### Define objective function
    # Flex Charging 
    #f_flex = cp.multiply(u , (TOU[:N_flex] - z[0]).reshape((N_flex))) # + cp.sum_squares(u) * lam_x
    f_flex = u.T @ TOU 
    g_flex = lam_h_c * 1 / z[2] 
    
    J_1 =  v[0] * (f_flex)
    
    # ASAP Charging
    f_asap = cp.sum(station_pow_max * (TOU[:N_asap] - z[1]))
    g_asap = lam_h_uc * 1 / z[2] 
    J_2 =  v[1] * (f_asap + g_asap)
    # Leave
    J_3 = cp.sum(TOU[:N_asap])

    J =   f_flex


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
    prob.solve(solver='GUROBI', warm_start=True)
    print("u:",np.round(u.value,1 ),"SOC:",np.round(SOC.value, 1))
    return  u.value, SOC.value

def plot_event(power_flex, SOC_flex, output_dir, filename, kwargs, ):

    plt.rcParams['figure.constrained_layout.use'] = True
    fig,(ax1, ax2, ax3, ax4) = plt.subplots(4,1,figsize=(12,9) ,)

    ts = np.arange(0,kwargs['N_flex'])
    # x_lim = max(ts)

    ax1.plot(ts, kwargs["TOU"])
    # ax1.set_xlim(x_lim)
    ax1.grid()
    ax1.set_xlabel("Interval")
    ax1.set_ylabel("TOU ($)")

    ax2.plot(ts, power_flex)
    # ax2.set_xlim(x_lim)
    ax2.grid()
    ax2.set_xlabel("Interval")
    ax2.set_ylabel("Power (kW)")
    ax2.set_title("Flex Charging Optimal Power")


    ax3.plot(ts[:kwargs["N_asap"]] , np.ones(shape=ts[:kwargs["N_asap"]].shape) * kwargs["station_pow_max"])
    # ax3.set_xlim(x_lim)
    ax3.grid()
    ax3.set_xlabel("Time")
    ax3.set_ylabel("Power (kW)")
    ax3.set_title("ASAP Charging Optimal Power")

    ax4.plot(SOC_flex)
    # ax4.set_xlim(x_lim)
    ax4.grid()
    ax4.set_xlabel("Time")
    ax4.set_ylabel("SOC (%)")
    ax4.set_title("Flex Charging SOC (%")
    plot_max_min_power([ax2,ax3], kwargs)

    plt.savefig(os.path.join(output_dir, filename))
    return 

def plot_max_min_power(ax_list, kwargs):
    for ax in ax_list:
        ax.plot(kwargs["station_pow_max"], label="Station Max Power (kW", color = 'red')
        # ax.plot(df_plot['Station Min Power (kW)'],label="Station Max Power (kW")))
    return 

def main():
    z = [50, 50, 10, 1]
    v = [0.45,0.45,0.1]


    kwargs = {
    "N_flex" : 40,
    "N_asap" : 10,
    "station_pow_max": 5,
    "lam_x" : 1,
    "lam_h_c" : 1,
    "lam_h_uc" : 1,
    "user_SOC_init"  :  0.3,
    "user_SOC_need": 0.8,
    "delta_k" : 1,
    "eff" : 1,
    "user_bat_cap" : 100}

    TOU = np.ones(shape = (kwargs["N_flex"],)) 
    TOU[10:15] = np.ones(shape = (5,)) * 5
    TOU[15:20] = np.ones(shape = (5,)) * 10
    TOU[20:30] = np.ones(shape = (10,)) * 5

    kwargs["TOU"] = TOU
    power_flex, SOC_flex = argmin_x(z,v, kwargs) 
    plot_event(power_flex, SOC_flex, "Figures/","test_3" , kwargs)

    return


if __name__ == "__main__":
    main()