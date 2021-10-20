import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from opt import Parameters, Problem, Optimization
import sys
import os
class Plot:
    def __init__(self, par):

        # TOU cost is hard coded in the Parameteters class with 30 min intervals 
        self.TOU_cost = par.TOU[0]
        self.intervals = np.arange( 0 , len(self.TOU_cost)) 

        # 

        ## 
        
    
    def event(self, output_dir = "Figures/",filename="fig.png"):
        """
        Function to plot single charging event: 
        1. TOU Cost 
        2. Power schedule for flex charging 
        3. Vehicle SOC, Flex 
        4. Power schedule for asap charging 
        5. Vehicle SOC, ASAP
        """
        ## Plot TOU 
        fig,(ax1) = plt.subplots(1,1,figsize=(10,3))
        ax1.plot(self.intervals, self.TOU_cost)
        ax1.set_xlim(self.intervals[0], self.intervals[-1])
        ax1.grid()
        ax1.set_xlabel("Interval")
        ax1.set_ylabel("TOU Cost ($")

        plt.savefig(os.path.join(output_dir,
        filename))

def main(new_event = False, time = None, pow_min = None, pow_max = None, overstay_duration = None, duration= None,batt_cap =None, SOC_need = None, SOC_init = None):
    par = Parameters()
    # if new_event:
    #     event = {"time" : time, "pow_min" : pow_min, "pow_max" : pow_max, "overstay_duration" : overstay_duration, "duration" : duration, "batt_cap" : batt_cap, "SOC_need" : SOC_need, "SOC_init" : SOC_init}
    #     prb = Problem(par = par, nargin = 1, event = event)
    # else:
    #     prb = Problem(par=par)
    # opt = Optimization(par, prb)
    # opt.run_opt()
    vis = Plot(par)
    
    vis.event(output_dir="Figures/", filename = "event_1.png")
    
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