import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from opt import Parameters, Problem, Optimization
import sys
import os
class Plot:
    
    def __init__(self, par, opt):

        # TOU cost is hard coded in the Parameteters class with 30 min intervals 
        self.num_interval= len(par.TOU[0])
        self.interval_size = par.Ts
        
        self.timerange = pd.date_range("1/1/2021", periods=self.num_interval, freq="{}H".format(self.interval_size ))
        self.TOU_cost =  pd.Series(par.TOU[0], index = self.timerange, name = "TOU Cost ($)")
        print(self.TOU_cost)

        # Event parameters: arrival_time, soc_init, soc_need, departure
        self.opt_time_start  = pd.Timestamp(year=2021, month=1, day=1, hour = int(opt.opt_time_start // 1), minute = int(60 *opt.opt_time_start % 1))
        self.opt_time_end_flex  =pd.Timestamp(year=2021, month=1, day=1, hour = int(opt.opt_time_end_flex  // 1),  minute = int(60 *opt.opt_time_end_flex % 1))
        self.opt_time_end_asap  = pd.Timestamp(year=2021, month=1, day=1, hour = int(opt.opt_time_end_asap // 1), minute = int(60 *opt.opt_time_end_asap % 1))

        self.flex_ts = pd.date_range(start=self.opt_time_start , end=self.opt_time_end_flex, freq="{}H".format(par.Ts))
        self.asap_ts = pd.date_range(start=self.opt_time_start , end=self.opt_time_end_asap, freq="{}H".format(par.Ts))
        
        # Optimization outputs 

        ## Question: Why flex output length mismatches the timelength? 
    
        self.opt_flex_powers = pd.Series(opt.opt_flex_powers.reshape(len(opt.opt_flex_powers),), index = flex_ts[2:], name = "Flew Power (kWh)")
        self.opt_asap_powers = pd.Series(opt.opt_asap_powers.reshape(len(opt.opt_asap_powers),), index = asap_ts, name = "ASAP Power (kWh)")
        self.opt_flex_SOCs = pd.Series(opt.pt_flex_SOCs.reshape(len(opt.opt_flex_SOCs),), index = flex_ts[2:], name = "Flex SOC (%)")

        # SOC? 
        print(self.opt_time_start)
        
    
    def event(self, output_dir = "Figures/",filename="fig.png", days = 2):
        """
        Function to plot single charging event: 
        2 Days 
        1. TOU Cost 
        2. Power schedule for flex charging 
        3. Vehicle SOC, Flex 
        4. Power schedule for asap charging 
        5. Vehicle SOC, ASAP

        What we want to check: 
        
        1. Is power schedule seem optimal? 
        2. max - min power limits 
        3. overstay? + charging and connection duration? 
        4. Vehicle SOC / Energy / Capacity? 
        """
        
        ts = pd.date_range(self.timerange[0], periods=self.num_interval * days, freq="{}H".format(self.interval_size ))
        
        df_plot = pd.DataFrame(index = ts, data = {"TOU Cost ($)": self.TOU_cost, "Flew Power (kWh)": self.opt_flex_powers, "ASAP Power (kWh)":self.opt_asap_powers})
        
        df_plot.to_csv("results.csv")
        print(self.TOU_cost)
        print(df_plot)
        ## Plot TOU 

        # fig,(ax1, ax2, ax3 ) = plt.subplots(3,1,figsize=(10,7))

        # ax1.plot()
        # ax1.set_xlim(self.timerange[0], self.timerange[-1])
        # ax1.grid()
        # ax1.set_xlabel("Time")
        # ax1.set_ylabel("TOU ($)")

        # ax2.plot(self.flex_ts[2:], self.opt_flex_powers)
        # ax2.set_xlim(self.timerange[0], self.timerange[-1])
        # ax2.set_xlabel("Time")
        # ax2.set_ylabel("Power (kWh)")



        # ax3.plot(self.asap_ts, self.opt_asap_powers)
        # ax3.set_xlim(self.timerange[0], self.timerange[-1])

        # plt.savefig(os.path.join(output_dir, filename))

def main(new_event = False, time = None, pow_min = None, pow_max = None, overstay_duration = None, duration= None,batt_cap =None, SOC_need = None, SOC_init = None):
    par = Parameters()
    if new_event:
        event = {"time" : time, "pow_min" : pow_min, "pow_max" : pow_max, "overstay_duration" : overstay_duration, "duration" : duration, "batt_cap" : batt_cap, "SOC_need" : SOC_need, "SOC_init" : SOC_init}
        prb = Problem(par = par, nargin = 1, event = event)
    else:
        prb = Problem(par=par)
    opt = Optimization(par, prb)
    opt.run_opt()
    vis = Plot(par,opt)
    
    # vis.event(output_dir="Figures/", filename = "event_1.png")
    
    return

if __name__ == "__main__":
    if len(sys.argv) > 2:
        time = float(sys.argv[1])
        pow_min = float(sys.argv[2])
        pow_max = float(sys.argv[3])
        overstay_duration = float(sys.argv[4])
        duration = float(sys.argv[5])
        batt_cap = float(sys.argv[6])
        SOC_need = float(sys.argv[7])
        SOC_init = float(sys.argv[8])
        main(new_event = True, time = time, pow_max = pow_max, pow_min = pow_min, overstay_duration= overstay_duration, duration= duration, batt_cap= batt_cap, SOC_init = SOC_init, SOC_need= SOC_need)
    else:
        print(sys.argv)
        main()