## The objective is to use CPLEX solver to solve single charger optimization 
import numpy as np
import datetime
from scipy import interpolate
import math
import timeit
import sys
# Need to install cplex
import cplex
from docplex.mp.model import Model
class Parameters:
    def __init__(self, sens_num_poles=np.arange(2, 18, 3), monte_num_sims=1, sim_starttime=7, sim_endtime=22,
                 sim_isFixedEventSequence=False, sim_isFixedSeed=False, sim_num_Events=50, Ts=0.25,
                 base_Tarriff_overstay=1.0, station_num_poles=8, eff=0.92, lam_x=10, lam_z_c=10, lam_z_uc=10,
                 lam_h_c=10, lam_h_uc=10, mu=1e4, soft_v_eta=1e-2, opt_eps=1e-4, VIS_DETAIL=True):
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
        # baseline  parameters
        self.base_tariff_overstay = base_Tarriff_overstay
        # TOU
        self.TOU = np.hstack((0.217 * np.ones((1, 34)),  # 0 - 8.5
                              0.244 * np.ones((1, 48 - 34)),  # 8.5 - 12
                              0.268 * np.ones((1, 72 - 48)),  # 12 - 16
                              0.244 * np.ones((1, 86 - 72)),  # 16 - 21.5
                              0.217 * np.ones((1, 96 - 86))))  # 22 - 24

        # charging station config
        self.station_num_poles = station_num_poles
        # number  of  charging  poles
        self.eff = eff  # power efficiency
        # dcm params
        self.dcm_choices = ['charging with flexibility', 'charging asap', 'leaving without charging']
        # pdfs
        self.pdf_visit = np.hstack((0.1 * np.ones((1, 7)),  # 0 - 7
                                    0.3 * np.ones((1, 5)),  # 7 - 12
                                    0.2 * np.ones((1, 2)),  # 12 - 14
                                    0.2 * np.ones((1, 2)),  # 14 - 16
                                    0.2 * np.ones((1, 6)),  # 16 - 22
                                    0.001 * np.ones((1, 2))))  # 22 - 24

        # regularization params
        self.lam_x = lam_x
        self.lam_z_c = lam_z_c
        self.lam_z_uc = lam_z_uc
        self.lam_h_c = lam_h_c
        # TODO: should be average overstay penalty in real data, should move to par
        self.lam_h_uc = lam_h_uc
        # TODO: should be average overstay penalty in real data, should move to par
        self.mu = mu
        self.soft_v_eta = soft_v_eta  # softening equality constraint for v; to avoid numerical error
        self.opt_eps = opt_eps
        # debug_mode
        self.VIS_DETAIL = VIS_DETAIL

class Problem:
    def __init__(self, nargin=0, par=Parameters(), **kwargs):
        self.par = par
        if nargin == 0:
            # user input
            self.user_time = 14.25
            self.user_SOC_init = 0.3
            self.user_SOC_need = 0.5
            self.user_batt_cap = 80  # kwh
            self.user_duration = 8  # hrs
            self.user_overstay_duration = 1
            self.station_pow_max = 7.2 # Station max power
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
        # dcm params
        # asc_flex = 2 + 0.2 * prb.user.duration;
        # asc_asap = 2.5;
        asc_flex = 2 + 0.401 * self.user_duration - 1.8531 * self.user_SOC_init  # 5.0583
        asc_asap = 1 + 0.865 * self.user_duration - 1.8531 * self.user_SOC_init  # 3.7088
        asc_leaving = 0
        energy_need = self.user_SOC_need * self.user_batt_cap
        # self.dcm_charging_flex.params = [-1 0 0 asc_flex].T
        # DCM parameters for choice 1 -- charging with flexibility
        # self.dcm_charging_asap.params = [0 - 1 0 asc_asap]';
        # DCM parameters for choice 2 -- charging as soon as possible

        self.dcm_charging_flex_params = np.array([[-0.1881 * energy_need], [0], [0], [asc_flex]])
        # % DCM parameters for choice 1 -- charging with flexibility
        self.dcm_charging_asap_params = np.array([[0], [- 0.1835 * energy_need], [0], [asc_asap]])
        # % DCM parameters for choice 2 -- charging as soon as possible
        self.dcm_leaving_params = np.array([[-0.01], [-0.01], [0.005], [asc_leaving]])
        # % DCM parameters for choice 3 -- leaving without charging
        self.THETA = np.vstack((self.dcm_charging_flex_params.T, self.dcm_charging_asap_params.T,
                                self.dcm_leaving_params.T))
        # problem specifications
        self.N_flex = int(self.user_duration / par.Ts)  # charging duration that is not charged, hour
        self.N_asap = math.floor((self.user_SOC_need - self.user_SOC_init) *
                                 self.user_batt_cap / self.station_pow_max / par.eff / par.Ts)
        self.TOU = interpolate.interp1d(np.arange(0, 24 - 0.25 + 0.1, 0.25), par.TOU, kind='nearest')(
            np.arange(self.user_time, 0.1 + self.user_time + self.user_duration - par.Ts, par.Ts)).T

class Optimization:
    def __init__(self, par, prb):
        self.Parameters = par
        self.Problem = prb
        self.opt_z = None
        self.opt_tariff_asap = None
        self.opt_tariff_flex = None
        self.opt_tariff_overstay = None

    def arg_v(self, z, x, **kwargs):
        """Choices"""
        """z: the price, x: state of charge"""
        """The Class Problem and Parameters bring values for this function to run properly"""
        conj = lambda v: np.sum([v[i]] * np.log(v[i]) for i in range(3))
        integ = lambda v: [v[i] * (self.Problem.THETA @ z)[i] for i in range(3)]
        obj = lambda v: ((np.sum(
            np.multiply(x[self.Problem.N_flex + 1:], (self.Problem.TOU[0:self.Problem.N_flex] - z[0])) + np.multiply(
                self.Parameters.lam_x, x[self.Problem.N_flex + 1:])) + self.Parameters.lam_z_c * (
                                      z[0] ** 2)) + self.Parameters.lam_h_c * 1 / z[2]) * v[0] + (np.sum(
            (self.Problem.station_pow_max * (self.Problem.TOU[0:self.Problem.N_asap] - z[1])) +
            self.Parameters.lam_z_uc * z[1] ** 2) + self.Parameters.lam_h_uc * 1 / z[2]) * v[1] + (
                            np.sum(self.Problem.station_pow_max * (self.Problem.TOU[0: self.Problem.N_asap] - 0))) * v[
                            2] + self.Parameters.mu * (conj(v) - np.sum(integ(v)))
        A = np.diag((-1 * np.ones((1, 3)))[0])
        b = np.zeros((3, 1))
        lower_bound = np.zeros((3, 1))
        beq = np.array([[-(1 - self.Parameters.soft_v_eta)], [1 + self.Parameters.soft_v_eta]])
        Aeq = np.vstack()
        model_v = Model("arg_v")
        ## Add all constraints
        model_v.add_constraint(b[0] - (A[0] @ b).item(0))
        model_v.add_constraint(b[1] - (A[1] @ b).item(0))
        model_v.add_constraint(b[2] - (A[2] @ b).item(0))
        model_v.add_constraint(beq[0] - (Aeq[0] @ b).item(0))
        model_v.add_constraint(beq[1] - (Aeq[1] @ b).item(0))
        model_v.add_constraint(lb=lower_bound[0])
        model_v.add_constraint(lb=lower_bound[1])
        model_v.add_constraint(lb=lower_bound[2])
        return model_v.minimize(obj)

    ##
    def arg_x(self, z, v):
        N = self.Problem.N_flex
        if sum(v) < 0 | (np.sum(v) < 1 - self.Parameters.soft_v_eta) | (np.sum(v) > 1 + self.Parameters.soft_v_eta):
            raise ValueError('[ ERROR] invalid $v$')
        model_x = Model("arg_x")
        ## Define the x variables
        x_vars = {(i, 1) for i in range(N + N + 1)}
        int1 = [x_vars[self.Problem.N_flex + 1:][i] * (self.Problem.TOU[:self.Problem.N_flex] - z[0])[i] for i in
                range(N)]
        int2 = [self.Parameters.lam_x * x_vars[self.Problem.N_flex + 1:][i] for i in range(N)]
        # The objective function
        obj_J = ((np.sum(int1) + np.sum(int2) + self.Parameters.lam_z_c * z[0] ** 2) + self.Parameters.lam_h_c * 1 / z[
            2]) * v[0] + (np.sum((self.Problem.station_pow_max * (
                    self.Problem.TOU[: self.Problem.N_asap] - z[1])) + self.Parameters.lam_z_uc * z[1] ** 2)
                          + self.Parameters.lam_h_uc * 1 / z[2]) * v[1] + np.sum(
            self.Problem.station_pow_max * (self.Problem.TOU[: self.Problem.N_asap] - 0)) * v[2]

        ## Prepare the data to properly put into the constraints
        A1L = np.hstack((np.zeros((1, N)), np.array([[-1]])))
        A1R = np.array(np.zeros((1, N)))
        A = np.hstack((A1L, A1R))
        b = -1 * self.Problem.user_SOC_need
        lb1 = np.zeros((N + 1, 1))
        lb1[-1] = self.Problem.user_SOC_need
        lb2 = self.Problem.station_pow_min * np.ones((N, 1))
        lb = np.vstack((lb1, lb2))
        ub1 = np.ones((N + 1, 1))  # state of charge
        ub2 = self.Problem.station_pow_max * np.ones((N, 1))  # max power
        ub = np.vstack((ub1, ub2))

        # equality constraints - system dynamics
        C1L = np.hstack((np.array([1]).reshape(1, 1), np.zeros((1, N))))  # initial soc for the first element
        C1R = np.zeros((1, N))

        C2L = np.hstack((np.diag((-1 * np.ones((1, N)))[0]), np.zeros((N, 1)))) + np.hstack(
            (np.zeros((N, 1)), np.diag(np.ones((1, N))[0])))
        C2R = np.diag((-1 * 1.47 * np.ones((1, N)))[0])

        d1 = self.Problem.user_SOC_init
        d2 = np.zeros((N, 1))
        Aeq = np.vstack((np.hstack((C1L, C1R)), np.hstack((C2L, C2R))))
        beq = np.vstack((d1, d2))
        model_x.add_constraint(A @ x_vars <= b)
        model_x.add_constraint(Aeq @ x_vars == beq)
        model_x.add_constraint(lb <= x_vars)
        model_x.add_constraint(x_vars <= ub)
        print(model_x)
        return model_x.minimize(obj_J)

    def arg_z(self, x, v):
        if sum(v) < 0 | (sum(v) < 1 - self.Parameters.soft_v_eta) | (sum(v) > 1 + self.Parameters.soft_v_eta):
            raise ValueError('[ ERROR] invalid {0}'.format(v))
        # Create a Model Z that will minimize the objective function
        model_z = Model("arg_z")
        # Put the data into proper format for the objective function
        int0 = lambda z: np.sum([self.Problem.THETA[0][i] * z[i] for i in range(4)])
        int1 = lambda z: np.sum([self.Problem.THETA[1][i] * z[i] for i in range(4)])
        int2 = lambda z: np.sum([self.Problem.THETA[2][i] * z[i] for i in range(4)])
        int3 = lambda z: np.sum([z[i] * (self.Problem.THETA.T @ v)[i] for i in range(4)])
        lse = lambda z: np.log(np.exp(int0(z)) + np.exp(int1(z)) + np.exp(int2(z)))
        J = lambda z: (((np.sum(
            (np.multiply(x[self.Problem.N_flex + 1:], (self.Problem.TOU[0: self.Problem.N_flex] - z[1]))) +
            np.multiply(self.Parameters.lam_x, x[self.Problem.N_flex + 1:])) + self.Parameters.lam_z_c * z[0] ** 2) +
                    self.Parameters.lam_h_c * 1 / z[2]) * v[0] + ((np.sum((self.Problem.station_pow_max * (
                    self.Problem.TOU[: self.Problem.N_asap] - z[1])) + self.Parameters.lam_z_uc * z[
                    1] ** 2) + self.Parameters.lam_h_uc * 1 / z[2]) * v[1])
                       + (np.sum(self.Problem.station_pow_max * (self.Problem.TOU[: self.Problem.N_asap])) * v[2])
                       + self.Parameters.mu * (lse(z) - int3(z)))
        A = [1, -1, 0, 0]
        Aeq = [0, 0, 0, 1]
        b = 0
        beq = 1
        lb = np.vstack((max(self.Problem.TOU) * np.ones((3,1)),0))
        ub = np.vstack((2*max(self.Problem.TOU) * np.ones((2,1)), 10 * max(self.Problem.TOU),1))
        # constraint1 = {'type': 'ineq', 'fun': lambda z: -(A[1] * z[1] + A[0] * z[0] + A[2] * z[2] + A[3] * z[3]) + b}

        constraint1 = -(A[1] + A[0] + A[2] + A[3]) + b

        # constraint2 = {'type': 'eq','fun': lambda z: Aeq[0] * z[0] + Aeq[1] * z[1] + Aeq[2] * z[2] + Aeq[3] * z[3] - beq}

        constraint2 = Aeq[0] + Aeq[1] + Aeq[2] + Aeq[3] - beq
        model_z.add_constraint(constraint1)
        model_z.add_constraint(constraint2)
        model_z.add_constraint(((lb[0], ub[0]), (lb[1], ub[1]), (lb[2], ub[2]), (lb[3], ub[3])))
        return model_z.minimize(J)

    def run_opt(self):
        start = timeit.timeit()
        def J_func(z, x, v):
            J = np.dot(np.array(
                [(np.sum((np.multiply(x[self.Problem.N_flex + 1:], (self.Problem.TOU[:self.Problem.N_flex] - z[0])))
                         + np.multiply(self.Parameters.lam_x, x[self.Problem.N_flex + 1:])) +
                  self.Parameters.lam_z_c * z[0] ** 2) + self.Parameters.lam_h_c * 1 / z[2],
                 np.sum((self.Problem.station_pow_max * (
                             self.Problem.TOU[0: self.Problem.N_asap] - z[1])) + self.Parameters.lam_z_uc * z[
                            1] ** 2) + self.Parameters.lam_h_uc * 1 / z[3],
                 np.sum(self.Problem.station_pow_max * (self.Problem.TOU[0: self.Problem.N_asap] - 0))]), v)
            return J

        itermax = 10000
        count = 0
        improve = np.inf
        zk = np.array([0, 0, 0, 1]).reshape(4, 1)
        # [z_c, z_uc, y, 1];
        xk = np.ones(
            (2 * self.Problem.N_flex + 1, 1))  # [soc0, ..., socN, u0, ..., uNm1]; - multiple dimensions 1 +  # of FLEX
        vk = np.array([0.45, 0.45, 0.1]).reshape(3, 1)  # [sm_c, sm_uc, sm_y]
        Jk = np.zeros((itermax, 1))
        while (count < itermax) & (improve >= 0) & (abs(improve) >= self.Parameters.opt_eps):
            count += 1
            z = zk
            x = xk
            v = vk
            Jk[count] = J_func(z, x, v)
            self.Problem.z0 = zk
            self.Problem.x0 = xk
            self.Problem.v0 = vk
            # update control variables
            zk = self.arg_z(xk, vk)
            xk = self.arg_x(zk, vk)
            vk = self.arg_v(zk, xk)

            # compute residual
            improve = Jk[count] - J_func(zk, xk, vk)

        self.opt_z = z
        self.opt_tariff_flex = z[0]
        self.opt_tariff_asap = z[1]
        self.opt_tariff_overstay = z[2]
        self.opt_x = x
        # update demand charge
        self.opt_peak_pow = max(x[self.Problem.N_flex + 2:])
        self.opt_flex_SOCs = x[0:self.Problem.N_flex + 1]
        self.opt_flex_powers = x[self.Problem.N_flex + 2:]
        self.opt_asap_powers = np.ones((self.Problem.N_asap, 1)) * self.Problem.station_pow_max
        self.opt_v = v
        self.opt_prob_flex = v[0]
        self.opt_prob_asap = v[1]
        self.opt_prob_leave = v[2]
        self.opt_J = Jk[0:count]
        self.opt_num_iter = count
        self.opt_prb = self.Problem
        self.opt_par = self.Parameters
        self.opt_time_start = self.Problem.user_time
        self.opt_time_end_flex = self.Problem.user_time + self.Problem.user_duration
        self.opt_time_end_asap = self.Problem.user_time + self.Problem.N_asap * self.Parameters.Ts
        end = timeit.timeit()
        if self.Parameters.VIS_DETAIL:
            print(
                ' {0} OPT DONE ({1} sec) sum(vk) = {2}, iterations = {3}\n'.format(datetime.datetime.now(), end - start,
                                                                                   sum(vk), count))

def main(new_event=False, time=None, pow_min=None, pow_max=None, overstay_duration=None, duration=None, batt_cap=None,
         SOC_need=None, SOC_init=None):
    par = Parameters()
    if new_event:
        event = {"time": time, "pow_min": pow_min, "pow_max": pow_max, "overstay_duration": overstay_duration,
                 "duration": duration, "batt_cap": batt_cap, "SOC_need": SOC_need, "SOC_init": SOC_init}
        prb = Problem(par=par, nargin=1, event=event)
    else:
        prb = Problem(par=par)
    opt = Optimization(par, prb)
    opt.run_opt()
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
        main(new_event=True, time=time, pow_max=pow_max, pow_min=pow_min, overstay_duration=overstay_duration,
             duration=duration, batt_cap=batt_cap, SOC_init=SOC_init, SOC_need=SOC_need)
    else:
        print(sys.argv)
        main()


