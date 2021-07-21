# Main credit for this Python script goes to Deep
import numpy as np
import cvxpy as cp
import datetime
from scipy import interpolate
import math
import timeit
from scipy.optimize import minimize
import sys
import pandas as pd
import matplotlib.pyplot as plt


class Parameters:
    def __init__(self,sens_num_poles = np.arange(2, 18, 3), monte_num_sims = 1, sim_starttime = 7, sim_endtime = 22,
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
        self.TOU = np.hstack((0.217 * np.ones((1, 34)), # 0 - 8.5
                   0.244 * np.ones((1, 48 - 34)), # 8.5 - 12
                   0.268 * np.ones((1, 72 - 48)), # 12 - 16
                   0.244 * np.ones((1, 86 - 72)), # 16 - 21.5
                   0.217 * np.ones((1, 96 - 86))))#22 - 24

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
    def __init__(self,nargin = 0, par = Parameters(),**kwargs):
        self.par = par
        if nargin == 0:
            # user input
            self.user_time = 14.25
            self.user_SOC_init = 0.3
            self.user_SOC_need = 0.5
            self.user_batt_cap = 80 # kwh
            self.user_duration = 8 # hrs
            self.user_overstay_duration = 1
            self.station_pow_max = 7.2
            self.station_pow_min = 0
        else:
            event = kwargs["event"]
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

    def argmin_v(self, z, x, **kwargs):
        lse_conj = lambda v: np.sum([v[i] * np.log(v[i]) for i in range(3)])
        int1 =  lambda v: [v[i] * (self.Problem.THETA @ z)[i] for i in range(3)]
        J = lambda v: ((np.sum(np.multiply(x[self.Problem.N_flex + 1 : ], (self.Problem.TOU[0:self.Problem.N_flex] - z[0])) + np.multiply(self.Parameters.lam_x, x[self.Problem.N_flex + 1 : ]))+  self.Parameters.lam_z_c * (z[0] ** 2)) + self.Parameters.lam_h_c * 1 / z[2]) * v[0] + (np.sum((self.Problem.station_pow_max * (self.Problem.TOU[0:self.Problem.N_asap] - z[1]))+
                    self.Parameters.lam_z_uc * z[1] ** 2) + self.Parameters.lam_h_uc * 1 / z[2]) * v[1] +  (np.sum(self.Problem.station_pow_max * (self.Problem.TOU[0: self.Problem.N_asap] - 0))) * v[2] + self.Parameters.mu * (lse_conj(v) - np.sum(int1(v)))
        A = np.diag((-1 * np.ones((1,3)))[0])
        b = np.zeros((3, 1))
        lb = np.zeros((3, 1))
        ub = np.array([1, 1, 0.3]).reshape(3,1)
        Aeq = np.vstack((-1 * np.ones((1, 3)), np.ones((1, 3))))
        beq = np.array([[-(1 - self.Parameters.soft_v_eta)], [1 + self.Parameters.soft_v_eta]])
        constraint1 = {'type': 'ineq', 'fun': lambda v:  b[0] - (A[0] @ v).item(0)}
        constraint2 = {'type': 'ineq', 'fun': lambda v:  b[1] - (A[1] @ v).item(0)}
        constraint3 = {'type': 'ineq', 'fun': lambda v:  b[2] - (A[2] @ v).item(0)}
        constraint4 = {'type': 'ineq', 'fun': lambda v:  beq[0] - (Aeq[0] @ v).item(0)}
        constraint5 = {'type': 'ineq', 'fun': lambda v:  beq[1] - (Aeq[1] @ v).item(0)}
        constraints = (constraint1, constraint2, constraint3, constraint4, constraint5)
        bounds = ((lb[0],ub[0]),(lb[1],ub[1]),(lb[2],ub[2]))
        v0 = self.Problem.v0
        res = minimize(J, v0, method='SLSQP', bounds=bounds, constraints=constraints)
        return res.x

    def argmin_x(self, z, v):
        N = self.Problem.N_flex

        if sum(v) < 0 | (np.sum(v) < 1 - self.Parameters.soft_v_eta) | (np.sum(v) > 1 + self.Parameters.soft_v_eta):
            raise ValueError('[ ERROR] invalid $v$')

        # cost    function
        # J = @(x) dot([sum((x(N + 2:end). * (prb.TOU(1:N) - z(1))). ^ 2) + par.lambda .h_c * 1 / z(3);
        # % sum((par.station.pow_max * (prb.TOU(1: N) - z(2))).^ 2) + par.lambda .h_uc * 1 / z(3);
        # % 1 / 3 * sum((par.station.pow_max * (prb.TOU(1: N) - 0)).^ 2)], v);
        x = cp.Variable(shape = (N + N + 1,1))
        int1 = [x[self.Problem.N_flex + 1 :][i] * (self.Problem.TOU[:self.Problem.N_flex] - z[0])[i] for i in range(N)]
        int2 = [self.Parameters.lam_x * x[self.Problem.N_flex + 1:][i] for i in range(N)]
        J = ((np.sum(int1) + np.sum(int2) + self.Parameters.lam_z_c * z[0] ** 2) + self.Parameters.lam_h_c * 1 / z[2]) * v[0] + (np.sum((self.Problem.station_pow_max * (self.Problem.TOU[: self.Problem.N_asap] - z[1])) + self.Parameters.lam_z_uc * z[1] ** 2)
           + self.Parameters.lam_h_uc * 1 / z[2]) * v[1] + np.sum(self.Problem.station_pow_max * (self.Problem.TOU[: self.Problem.N_asap ] - 0)) * v[2]

        # inequality constraint
        A1L = np.hstack((np.zeros((1, N)), np.array([[-1]])))
        A1R = np.array(np.zeros((1, N)))
        A = np.hstack((A1L , A1R))
        b = -1 * self.Problem.user_SOC_need

        # lower bound - power min
        lb1 = np.zeros((N + 1, 1))
        lb1[-1] = self.Problem.user_SOC_need
        lb2 = self.Problem.station_pow_min * np.ones((N,1))
        lb = np.vstack((lb1, lb2))

        # upper bound - power max
        ub1 = np.ones((N + 1, 1)) # soc
        ub2 = self.Problem.station_pow_max * np.ones((N,1)) # power
        ub = np.vstack((ub1, ub2))

        # equality constraints - system dynamics
        C1L = np.hstack((np.array([1]).reshape(1,1) ,np.zeros((1, N)))) # initial soc for the first element
        C1R = np.zeros((1, N))

        C2L = np.hstack((np.diag((-1 * np.ones((1, N)))[0]),np.zeros((N, 1)))) + np.hstack((np.zeros((N, 1)), np.diag(np.ones((1, N))[0])))
        C2R = np.diag((-1 * self.Parameters.eff* self.Parameters.Ts/ self.Problem.user_batt_cap * np.ones((1, N)))[0])

        d1 = self.Problem.user_SOC_init
        d2 = np.zeros((N, 1))

        Aeq = np.vstack((np.hstack((C1L, C1R)), np.hstack((C2L, C2R))))
        beq = np.vstack((d1, d2))
        constraints = [A @ x <= b, Aeq @ x == beq, lb <= x, x <= ub]
        obj = cp.Minimize(J)
        x.value = self.Problem.x0
        problem = cp.Problem(obj, constraints)
        problem.solve(solver="ECOS", warm_start=True)

        return x.value

    def argmin_z(self, x, v):

        if sum(v) < 0 | (sum(v) < 1 - self.Parameters.soft_v_eta) | (sum(v) > 1 + self.Parameters.soft_v_eta):
            raise ValueError('[ ERROR] invalid {0}'.format(v))
        int0 = lambda z: np.sum([self.Problem.THETA[0][i] * z[i] for i in range(4)])
        int1 = lambda z: np.sum([self.Problem.THETA[1][i] * z[i] for i in range(4)])
        int2 = lambda z: np.sum([self.Problem.THETA[2][i] * z[i] for i in range(4)])
        lse = lambda z: np.log(np.exp(int0(z)) + np.exp(int1(z)) + np.exp(int2(z)))
        int3 = lambda z: np.sum([z[i] * (self.Problem.THETA.T @ v)[i] for i in range(4)])
        J = lambda z: (((np.sum((np.multiply(x[self.Problem.N_flex + 1 :], (self.Problem.TOU[0: self.Problem.N_flex] - z[0]))) +
                     np.multiply(self.Parameters.lam_x, x[self.Problem.N_flex + 1:])) +
              self.Parameters.lam_z_c * z[0] ** 2) +
             self.Parameters.lam_h_c * 1 / z[2]) * v[0] + ((np.sum((self.Problem.station_pow_max * (self.Problem.TOU[: self.Problem.N_asap ] - z[1]))+self.Parameters.lam_z_uc * z[1] ** 2) + self.Parameters.lam_h_uc * 1 / z[2]) * v[1])
        + (np.sum(self.Problem.station_pow_max * (self.Problem.TOU[: self.Problem.N_asap])) * v[2]) \
        + self.Parameters.mu * (lse(z) - int3(z)))
        A = [1, -1, 0, 0]
        b = 0
        #charging tariff for charging asap must be bigger

        # lower and upper bound
        lb = np.vstack((max(self.Problem.TOU) * np.ones((3, 1)), 0))
        ub = np.vstack((2 * max(self.Problem.TOU) * np.ones((2, 1)), 10 * max(self.Problem.TOU), 1))
        #equality constraint
        Aeq = [0, 0, 0, 1]
        beq = 1
        constraint1 = {'type': 'ineq', 'fun': lambda z: -(A[1] * z[1] + A[0] * z[0] + A[2] * z[2] + A[3] * z[3]) + b}
        constraint2 = {'type': 'eq', 'fun': lambda z: Aeq[0] * z[0] + Aeq[1] * z[1] + Aeq[2] * z[2] + Aeq[3] * z[3] - beq}
        constraints = (constraint1, constraint2)
        bounds = ((lb[0], ub[0]), (lb[1], ub[1]), (lb[2], ub[2]), (lb[3], ub[3]))
        z0 = self.Problem.z0
        res = minimize(J, z0, method='SLSQP', bounds=bounds, constraints=constraints)
        return res.x

    def run_opt(self):
        start = timeit.timeit()
        def J_func(z, x, v):
            J = np.dot(np.array([(np.sum((np.multiply(x[self.Problem.N_flex + 1:], (self.Problem.TOU[:self.Problem.N_flex ] - z[0])))
                                         + np.multiply(self.Parameters.lam_x, x[self.Problem.N_flex + 1:]))+
                             self.Parameters.lam_z_c * z[0] ** 2) + self.Parameters.lam_h_c * 1 / z[3],
                         np.sum((self.Problem.station_pow_max * (self.Problem.TOU[0: self.Problem.N_asap ] - z[1])) + self.Parameters.lam_z_uc * z[1] ** 2) + self.Parameters.lam_h_uc * 1 / z[3],
                        np.sum(self.Problem.station_pow_max * (self.Problem.TOU[0: self.Problem.N_asap ] - 0))]), v)
            return J
        itermax = 10000
        count = 0
        improve = np.inf
        zk = np.array([0, 0, 0, 1]).reshape(4,1)
        # [z_c, z_uc, y, 1];
        xk = np.ones((2 * self.Problem.N_flex + 1, 1)) # [soc0, ..., socN, u0, ..., uNm1]; - multiple dimensions 1 +  # of FLEX
        vk = np.array([0.45, 0.45, 0.1]).reshape(3,1)                     # [sm_c, sm_uc, sm_y]
        Jk = np.zeros((itermax, 1))
        while (count < itermax) & (improve >= 0) & (abs(improve) >= self.Parameters.opt_eps):
            z = zk
            x = xk
            v = vk
            Jk[count] = J_func(z, x, v)
            self.Problem.z0 = zk
            self.Problem.x0 = xk
            self.Problem.v0 = vk
            #update control variables
            zk = self.argmin_z(xk, vk)
            xk = self.argmin_x(zk, vk)
            vk = self.argmin_v(zk, xk)

            # compute residual
            improve = Jk[count] - J_func(zk, xk, vk)
            count += 1
            print(count)
        opt = {}
        opt["z"] = zk
        opt["tariff_flex"] = zk[0]
        opt["tariff_asap"] = zk[1]
        opt["tariff_overstay"] = zk[2]
        opt["x"] = xk
        # update demand charge
        opt["peak_pow"] = max(xk[self.Problem.N_flex + 1:])
        opt["flex_SOCs"] = xk[0:self.Problem.N_flex + 1]
        opt["flex_powers"] = xk[self.Problem.N_flex + 1:]
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
        if self.Parameters.VIS_DETAIL:
            print(
                ' {0} OPT DONE ({1} sec) sum(vk) = {2}, iterations = {3}\n'.format(datetime.datetime.now(), end - start,
                                                                                   sum(vk), count))
        return opt

class Simulation:
    def __init__(self, par = Parameters(), t = None):
        self.par = par
        t = t if t is not None else np.arange(par.sim_starttime, par.sim_endtime + 0.1, par.Ts)
        self.t = t
        self.power = np.zeros(t.shape)
        self.profit = np.zeros(t.shape)
        self.profit_charging_uc = np.zeros(t.shape)
        self.profit_charging_c = np.zeros(t.shape)
        self.profit_overstay = np.zeros(t.shape)
        self.occ_total = np.zeros(t.shape)
        self.occ_empty = np.zeros(t.shape)
        self.occ_charging = np.zeros(t.shape)
        self.occ_overstay = np.zeros(t.shape)
        self.overstay_duration = np.zeros(t.shape)
        choice = np.zeros((self.par.sim_num_events, 1))
        self.choice = choice
        choice_probs = np.zeros((self.par.sim_num_events, 3))
        self.choice_probs = choice_probs
        self.num_service = np.zeros(t.shape)
        self.tot_decision = 0
        self.t = t
        self.control = np.zeros((len(t), 3))
        self.opts = {}

    def gen_events_one_day(self, param= None, seed = None):
        if param is None:
            par = self.par
        elif seed is None:
            par = param
        else:
            par = param
            seed_val = seed
            np.random.seed(seed_val)
        if par.sim_isFixedEventSequence:
            act_data = pd.read_csv("real_act_data_1day.csv")
            num_events = len(act_data)
            event_idx = range(num_events)
        else:
            num_events = par.sim_num_events
            act_data = pd.read_csv("real_act_data.csv")
            event_idx = np.sort(np.random.randint(0, len(act_data) - 1, num_events))
        self.events = {}
        self.events["inp"] = []
        self.events["time"] = []
        self.events["triggered"] = np.full((num_events, 1), False)

        for i in range(num_events):
            n = event_idx[i]
            if act_data.loc[n, "duration"] < 0.3:
                continue
            event = {}
            event["time"] = act_data.loc[n, "arr_time"]
            event["SOC_init"] = act_data.loc[n, "soc_init"]
            event["SOC_need"] = act_data.loc[n, "soc_need"]
            event["batt_cap"] = act_data.loc[n, "batt_cap"]
            event["duration"] = act_data.loc[n, "duration"]
            event["overstay_duration"] = act_data.loc[n, "overstay_duration"]
            event["pow_max"] = act_data.loc[n, "power_max"]
            event["pow_min"] = 0
            self.events["inp"].append(event)
            self.events["time"].append(event["time"])
        return

    def run_sim_one_day(self, params = None, events = None, vis = False):
        if params is None:
            par = self.par
            self.gen_events_one_day(param = par)
        elif events is None:
            par = params
            self.gen_events_one_day(param = par)
        else:
            par = params
            self.events = events

        if par.VIS_DETAIL:
            print("[{0} INIT] INITIALIZATION DONE\n".format(datetime.datetime.now()))

        t = np.arange(par.sim_starttime, par.sim_endtime + 0.1, par.Ts)
        i_k = 0
        i_event = 0
        station = {}
        station["num_occupied_pole"] = 0
        station["num_empty_pole"] = par.station_num_poles
        for k in t:
            if i_event < len(self.events["time"]):
                if np.any([round(i/par.Ts)*par.Ts == k for i in self.events["time"]]):
                    inds_events_k = [v for v, j in enumerate(self.events["time"]) if round(j/par.Ts)*par.Ts == k]
                    for j in range(len(inds_events_k)):
                        if (self.events["inp"][inds_events_k[j]]["duration"] <= par.sim_endtime - k) & (station["num_empty_pole"] > 0):
                            self.tot_decision  += 1
                            self.events["triggered"][i_event] = True
                            prb = Problem(event = self.events["inp"][inds_events_k[j]])
                            optim = Optimization(prb= prb, par=par)
                            opt = optim.run_opt()
                            self.opts[i_event] = opt
                            rc = np.random.uniform(0, 1)
                            if rc <= opt["prob_flex"]:
                                opt["choice"] = 0
                                opt["time_end"] = opt["time_end_flex"]
                                opt["powers"] = opt["flex_powers"]
                                opt["price"] = opt["tariff_flex"]
                            elif rc <= opt["prob_flex"] + opt["prob_asap"]:
                                opt["choice"] = 1
                                opt["time_end"] = opt["time_end_asap"]
                                opt["powers"] = opt["asap_powers"]
                                opt["price"] = opt["tariff_asap"]
                            else:
                                opt["choice"] = 2
                            self.choice_probs[i_event][:] = opt["v"]
                            self.choice[i_event] = opt["choice"]
                            self.control = opt["z"][0:2]
                            if par.VIS_DETAIL:
                                print('[{0} EVENT] time = {1}, CHOICE = {2}\n'.format(datetime.datetime.now(), k, par.dcm_choices[opt["choice"]]))

                            if opt["choice"] <= 1:
                                opt["time_leave"], duration = self.get_rand_os_duration(optim = opt, prb = prb)
                                self.overstay_duration[i_k] += duration
                                self.num_service[i_k] += 1
                                station["num_occupied_pole"] += 1
                                station["num_empty_pole"] -= 1
                                station["EV" + str(self.tot_decision)] = opt
                            else:
                                if par.VIS_DETAIL:
                                    if station["num_empty_pole"] == 0:
                                        print("[{0} EVENT] SKIPPED (event {1}) due to full occupancy\n".format(datetime.datetime.now(),i_event))
                                    else:
                                        print("[{0} EVENT] SKIPPED (event {1}) due to violating operating hours\n".format(datetime.datetime.now(), i_event))
                        i_event = i_event + 1
            if len(station.keys()) > 0:
                keys_to_delete = []
                for ev in station.keys():
                    pop = False
                    if "EV" in ev:
                        if k < station[ev]["time_end"]:
                            TOU = interpolate.interp1d(np.arange(0, 24 - 0.25 + 0.1, 0.25), par.TOU, kind = 'nearest')(k).T
                            opt = station[ev]
                            dur = np.arange(opt["time_start"],opt["time_end"],par.Ts)
                            if len(dur) > 1  & len(opt["powers"]) > 1:
                                power = interpolate.interp1d(dur, opt["powers"])(k).T
                            else:
                                power = opt["powers"][0]

                            if opt["choice"] == 1:
                                self.profit_charging_uc[i_k] += par.Ts * power * (station[ev]["price"] - TOU)
                            else:
                                self.profit_charging_c[i_k] += par.Ts * power * (station[ev]["price"] - TOU)

                            self.power[i_k] += power
                            self.occ_charging[i_k] += 1
                        else:
                            if k < station[ev]["time_end"]:
                                self.profit_overstay[i_k] += par.Ts * station[ev]["tariff_overstay"]
                                self.occ_overstay[i_k] += 1
                            else:
                                pop = True
                                keys_to_delete.append(ev)
                                print("POP HAPPENED")
                                station["num_occupied_pole"] -= 1
                                station["num_empty_pole"] += 1
                    elif "occ" in ev:
                        self.occ_total[i_k] = station["num_occupied_pole"]
                        self.occ_empty[i_k] = station["num_empty_pole"]
                if len(keys_to_delete) > 0:
                    for i in keys_to_delete:
                        del station[i]
            i_k += 1
        self.par = par
        if vis:
            self.vis_sim_one_day(output = False)

    def get_rand_os_duration(self, optim = None, prb = Problem()):
        scale = 100
        lam = prb.user_overstay_duration * scale * (self.par.base_tariff_overstay / optim["tariff_overstay"])
        ran = np.arange(0, scale * 10 + 1)
        fact = [math.factorial(i) for i in ran]
        pdf = []
        for i in range(len(ran)):
            try:
                pdf.append(np.exp(-lam) * np.power(lam, ran[i]) / fact[i])
            except:
                pdf.append(np.nan)
        cdf = np.cumsum(pdf)
        r = (1 - min(cdf)) * np.random.random_sample() + min(cdf)
        duration = 0
        try:
            duration = interpolate.interp1d(cdf, ran)(r)/scale
        except:
            duration = ran[[i for i,j in enumerate(cdf) if j >= r][0]]/scale
        # if duration.shape == (1,):
        #     duration = 32
        overstay_endtime = optim["time_end"] + duration
        if np.isnan(duration):
            ValueError("Error: NAN Duration")
        return overstay_endtime, duration


    def vis_sim_one_day(self, output = False, output_dict = None):
        options = {"display": True, "temporals": True, "choices": True}
        if output:
            return options

        if output_dict is not None:
            options["display"] = output_dict["display"]
            options["temporals"] = output_dict["temporals"]
            options["choices"] = output_dict["choices"]

        if options["display"]:
            if options["temporals"]:
                fig,(ax1,ax2,ax3,ax4,ax5) = plt.subplots(5,1,figsize=(17,9))
                ax1.plot(self.t, self.power, linewidth = 1.5)
                ax1.set_xlim(self.t[0], self.t[-1])
                ax1.grid()
                ax1.set_xlabel("Hour of the Day")
                ax1.set_ylabel("Power (kW)")

                plt.rcParams.update({'font.size': 15, "legend.fontsize": 15, 'legend.handlelength': 2})

                ax2.plot(self.t,self.profit_charging_c + self.profit_charging_uc + self.profit_overstay, linewidth = 1.5, label = "Instant")
                tot_profit = np.cumsum(self.profit_charging_uc + self.profit_charging_c + self.profit_overstay)
                ax2.plot(self.t, tot_profit, linewidth = 1.5, label = "Net")
                ax2.set_xlim(self.t[0],self.t[-1])
                ax2.grid(True)
                ax2.set_xlabel("Hour of the Day")
                ax2.set_ylabel("Profit ($)")
                ax2.legend()

                ax3.plot(self.t, self.occ_total, linewidth = 1.5, label = "Total")
                ax3.plot(self.t, self.occ_overstay, linewidth = 1.5, label = "Overstay")
                ax3.plot(self.t, self.occ_charging, linewidth = 1.5, label = "Charging")
                ax1.set_xlim(self.t[0],self.t[-1])
                ax3.legend()
                ax3.grid()
                ax3.set_xlabel("Hour of the day")
                ax3.set_ylabel("# of Vehicles")

                ax4.plot(self.t, np.cumsum(self.overstay_duration), linewidth = 1.5)
                ax4.set_xlim(self.t[0], self.t[-1])
                ax4.grid(True)
                ax4.set_xlabel("Hour of the day")
                ax4.set_ylabel("Overstay Duration (hours)")

                ax5.stem(self.t, self.num_service, label = "Instant")
                ax5.plot(self.t, np.cumsum(self.num_service), linewidth = 1.5, label = "Net")
                ax5.set_xlim(self.t[0], self.t[-1])
                ax5.legend()
                ax5.set_xlabel("Hour of the day")
                ax5.grid(True)
                ax5.set_ylabel("# of service")
                plt.savefig("one_day_temporal.eps")
                plt.savefig("one_day_temporal.png")

            if options["choices"]:
                fig = plt.figure(figsize = (10,5))
                choices = [j for i,j in enumerate(self.choice) if ~np.isnan(self.choice)[i]]
                choice_probs = [self.choice_probs[i] for i in range(len(self.choice_probs)) if any(~np.isnan(self.choice_probs)[i])]
                choice_labels = np.zeros((len(choices),1))
                choice_times = [j for i,j in enumerate(self.events["time"]) if self.events["triggered"][i] == True]
                choice_times_str = [str(round(choice_times[i] * 100)/100) for i in range(len(choice_times))]
                choice_times_str_filtered = []
                for i in range(len(choice_times_str)):
                    choice_times_str_filtered.append(['(' + choice_times_str[i].strip() + ')'])
    
                for j in range(len(choice_probs)):
                    choice_labels[j] = [1 if choices[j] == 0 else 0][0] * 1/2 * choice_probs[j][0] +\
                                       [1 if choices[j] == 1 else 0][0] * (choice_probs[j][0]
                                                                           + 0.5 * choice_probs[j][1]) +\
                                       [1 if choices[j] == 2 else 0][0] * (choice_probs[j][0] + choice_probs[j][1] + 0.5 * choice_probs[j][2])
                for i in range(3):
                    plt.fill_between(range(len(choice_labels)),[choice_probs[j][i] for j in range(len(choice_probs))])
                plt.annotate("charging_flex",(1, 0.04))
                plt.annotate("charging asap",(1, sum(choice_probs[0][:2]) - 0.04))
                plt.annotate("leaving without charging", (.1, sum(choice_probs[0][:3]) - 0.04))
                plt.scatter(range(len(choice_labels)),choice_labels, color = "black")
                plt.xlabel("event #")
                plt.ylabel("Probability [0 1]")
                plt.rc('figure', titlesize= 1)
                plt.title("Choice Probabilities and Decisions \n total: {0}, flex: {1}, asap: {2}, leave: {3}".format(len(choices), sum([1 for i in range(len(choices)) if choices[i]==0]), sum([1 for i in range(len(choices)) if choices[i]==1]), sum([1 for i in range(len(choices)) if choices[i]==2])))
                plt.xlim(1, len(choice_probs))
                plt.ylim(0,1)
                plt.legend("", "choice")
                plt.savefig("one_day_choice.png")
                plt.savefig("one_day_choice.eps")

                # fig, (ax2,ax3) = plt.subplots(1,2)
                # ax1.fill_between(range(len(choice_labels)), choice_probs)
                # ax1.annotate("charging_flex", 1.5, 0.04)
                # ax1.annotate("charging asap", 1.5, sum(choice_probs[0][:2]) - 0.04)
                # ax1.annotate("leaving without charging", 1.5, sum(choice_probs[0][:3]) - 0.04)
                # ax1.scatter(range(len(choice_labels)), choice_labels, color="black")
                # ax1.xlabel("event #")
                # ax1.ylabel("Probability [0 1]")
                # ax1.title("Choice Probabilities and Decisions \n total: {0}, flex: {1}, asap: {2}, leave: {3}".format(
                #     len(choices), sum(choices == 0), sum(choices == 1), sum(choices == 2)))
                # ax1.xlim(1, len(choice_probs))
                # ax1.ylim(0, 1)


def main(new_event = False, time = None, pow_min = None, pow_max = None, overstay_duration = None, duration= None,batt_cap =None, SOC_need = None, SOC_init = None):
    par = Parameters()
    prob = Problem()
    sim = Simulation(par = par)
    sim.run_sim_one_day(vis=True)
    return


if __name__ == "__main__":
    if len(sys.argv) > 2:
        time = float(sys.argv[1])
        print(time)
        print(type(time))
        print(sys.argv[1])
        print(type(sys.argv[1]))
        pow_min = float(sys.argv[8])
        pow_max = float(sys.argv[7])
        overstay_duration = float(sys.argv[6])
        duration = float(sys.argv[5])
        batt_cap = float(sys.argv[4])
        SOC_need = float(sys.argv[3])
        SOC_init = float(sys.argv[2])
        main(new_event = True, time = time, pow_max = pow_max, pow_min = pow_min, overstay_duration= overstay_duration, duration= duration, batt_cap= batt_cap, SOC_init = SOC_init, SOC_need= SOC_need)
    else:
        print(sys.argv)
        main()
