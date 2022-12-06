# PEV-OPT-EXPERIMENT 
Python implementation of the joint price and power optimizer and analysis of the experimental results. 

Optimizer: Python_Code/optimizer.py

Demo: DEMO.ipynb

KPI Analysis: kpi_analysis.ipnyb



## Workflow

![demo](demo/decision_flow.png#style=centerme) 
<center><h5><em>Fig. 1. PEV charging station work flow: the decision process when a user plug-in a PEV</em></h5></center>

Fig. 1 illustrates the charger operation for a single PEV driver (denoted as “user”). Upon arrival to the PEV charging station, the user inputs the following information: 
- intended parking duration
- desired added range in miles.

In sequence, the user receives the pricing for two charging service options in \[$/kW], and an overstay penalty in \[$/hour] (bottom box in Fig. 1). These prices are computed by the pricing policy controller. Given the prices, the user chooses one of the following three options:
- charging-flexibility (controlled charging, time flexibility granted by customer): The needed energy is guaranteed to be delivered upon departure. However, the station operator may optimize the charging schedule.
- charging-asap (uncontrolled charging, no time flexibility permitted): The PEV is charged at the max power continuously, starting immediately, until the vehicle departs or the battery is full.
- leave: the driver leaves the station without charging. 

If a charger is vacant and the user decides for either charging  service (charging-flexibility or charging-asap), then  the charger will be occupied for the entire parking duration. When the user departs, the user pays the service fee (including overstay fees if applicable). The charger then becomes available to others. If the user decides to leave withoutcharging (leave), the charger remains open to others.

## Controller 
There are two strategies of optimization: (i) single-charger optimization \[1] and (ii) station-wide optimization \[2]. The single-charger optimization strategy only considers one charger (which the user is currently occupying) in optimal charging schedule. The station-wide optimization strategy considers all chargers (which the previous users are currently occupying) in optimal charging schedule. 



## References
\[1] T. Zeng\*, S. Bae\*, B. Travacca, S.J. Moura, _"Inducing Human Behavior to Maximize Operation Performance at PEV Charging Station"_, IEEE Trans. on Smart Grid. Availabel here: https://doi.org/10.1109/TSG.2021.3066998. \*equal

\[2] S. Bae, T. Zeng, B.  Travacca, and S. Moura, _Inducing Human Behavior to Alleviate Overstay at PEV Charging Station_, American Control Conference, 2020. [preprint](https://arxiv.org/pdf/1912.02341.pdf).
```
@article{zeng2021inducing,
  title={Inducing Human Behavior to Maximize Operation Performance at PEV Charging Station},
  author={Zeng, Teng and Bae, Sangjae and Travacca, Bertrand and Moura, Scott},
  journal={IEEE Transactions on Smart Grid},
  year={2021},
  publisher={IEEE}
}

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
* Ayse Tugba Ozturk: tugbaozturk@berkeley.edu
* Sangjae Bae: sangjae.bae@berkeley.edu
* Teng Zeng: tengzeng@berkeley.edu
