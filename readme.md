# PEV-OPERATION-HIL
Plug-in Electric Vehicle Charging Station Operation Simulator -- Human-In-the-Loop. The simulator implements the online controller that finds optimal charging and pricing policy. 

## Workflow
Fig.  1  illustrates  the  charger  operation  for  a  single  PEV driver (denoted as “user”). Upon arrival to the PEV charging station, the user inputs the following information: 
- intended parking duration
- desired added range in miles.

In sequence, the user receives the pricing for two charging service options in \[$/kW], and an overstay penalty in \[$/hour] (bottom box in Fig. 1). These prices are computed by the pricing policy controller. Given the prices, the user chooses one of the following three options:
- charging-flexibility (controlled charging, time flexibility granted by customer): The needed energy is guaranteed to be delivered upon departure. However, the station operator may optimize the charging schedule.
- charging-asap (uncontrolled charging, no time flexibility permitted): The PEV is charged at the max power continuously, starting immediately, until the vehicle departs or the battery is full.
- leave: the driver leaves the station without charging. 

If a charger is vacant and the user decides for either charging  service (charging-flexibility or charging-asap), then  the charger will be occupied for the entire parking duration. When the user departs, the user pays the service fee (including overstay fees if applicable). The charger then becomes available to others. If the user decides to leave withoutcharging (leave), the charger remains open to others.

![demo](demo/decision_flow.png#style=centerme)
<center>Fig. 1. PEV charging station work flow: the decision process when a user plug-in a PEV</center>

## Controller 
There are two strategies of optimization: (i) single-charger optimization \[1] and (ii) station-wide optimization \[2]. The single-charger optimization strategy only considers one charger (which the user is currently occupying) in optimal charging schedule. The station-wide optimization strategy considers all chargers (which the previous users are currently occupying) in optimal charging schedule. 

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


## References
\[1] S. Bae, T. Zeng, B.  Travacca, and S. Moura, _Inducing Human Behavior to Alleviate Overstay at PEV Charging Station_, American Control Conference, 2020. [preprint](https://arxiv.org/pdf/1912.02341.pdf).

\[1] T. Zeng\*, S. Bae\*, B. Travacca, S.J. Moura, _"Inducing Human Behavior to Maximize Operation Performance at PEV Charging Station"_, In preparation.\*equal
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

## Contact
* Sangjae Bae: sangjae.bae@berkeley.edu
* Teng Zeng: hustlejanton@berkeley.edu
