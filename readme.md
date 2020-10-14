# PEV-OPERATION-HIL
Plug-in Electric Vehicle Charging Station Operation Simulator -- Human-In-the-Loop. The simulator implements the online controller that finds optimal charging and pricing policy. 

![demo](demo/decision_flow.png)

* S. Bae, T. Zeng, B.  Travacca, and S. Moura, _Inducing Human Behavior to Alleviate Overstay at PEV Charging Station_, American Control Conference, 2020. [preprint](https://arxiv.org/pdf/1912.02341.pdf).

* T. Zeng\*, S. Bae\*, B. Travacca, S.J. Moura, _"Inducing Human Behavior to Maximize Operation Performance at PEV Charging Station"_, In preparation.\*equal
```
@inproceedings{bae2020inducing,
  title={Inducing Human Behavior to Alleviate Overstay at PEV Charging Station},
  author={Bae, Sangjae and Zeng, Teng and Travacca, Bertrand and Moura, Scott},
  booktitle={2020 American Control Conference (ACC)},
  pages={2388--2394},
  year={2020},
  organization={IEEE}
}
```


## Prerequisites
- Matlab R2014a or higher
- Parallel Computing Toolbox (for sensitivity analysis)


## Setup
In the Matlab command window
```bash
cd ~/pev-operation-hil
run setup
```

## Project structure
- ```scripts```: contains matlab script files
    - ```init_par``` : initialize par structure that stores constant (hyper) parameters -- **configuration** in this file    
    - ```init_prb``` : initialize prb structure that stores event specific parameteres 
    - ```run_sim_one_day.m``` : run simulation of one day operation
    - ```run_sim_one_event.m``` : run simulation of one control event
    - ```run_sim_monte.m``` : run monte carlo simulations
    - ```run_sensitivity_analysis.m``` : run sensitivity analysis (by the number of charging poles)
    - ```get_rand_os_duration.m``` : randomly sample an overstay duration
    - ```init_sim``` : initialize sim structure that stores one day operation simulation result
    - ```run_opt``` : solve optimization with block coordinate descent (BCD) algorithm for a current EV
    - ```argmin_z``` : update zk in BCD algorithm
    - ```argmin_x``` : update xk in BCD algorithm
    - ```argmin_v``` : update vk in BCD algorithm
    - ```constr_J``` : cost function for the station-wide optimization
    - ```vis_sim_monte.m``` : visualize a result of monte carlo simulations results. Run this file to open a file brower and select a .mat file in /monte-sim-results 
    - ```vis_sim_one_day.m``` : visualize a simulation result of one day operation
    - ```vis_sim_one_event.m``` : visualize a simulation result of one control event
    - ```*_station.m``` : extension with ```_station``` indicates the station-wide optimization
