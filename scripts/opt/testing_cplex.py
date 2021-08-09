## This is not fully complete since we need to run the main optimization function 

import numpy as np
import datetime
from scipy import interpolate
import math
import timeit
from scipy.optimize import minimize
import sys
import cplex 
from docplex.mp.model import Model 

class Optimization:
    def __init__(self, par, prb):
        self.Parameters = par
        self.Problem = prb
        self.opt_z = None
        self.opt_tariff_asap = None
        self.opt_tariff_flex = None
        self.opt_tariff_overstay = None
    def arg_v(self, z, x, **kwargs):
        """z: the price, x: state of charge"""
        """The Class Problem and Parameters bring values for this function to run properly"""
        conj = lambda v: np.sum([v[i]] * np.log(v[i]) for i in range(3))
        integ = lambda v: [v[i]*(self.Problem.THETA @ z)[i] for i in range(3)]
        obj = lambda v: ((np.sum(np.multiply(x[self.Problem.N_flex + 1 : ], (self.Problem.TOU[0:self.Problem.N_flex] - z[0])) + np.multiply(self.Parameters.lam_x, x[self.Problem.N_flex + 1 : ]))+  self.Parameters.lam_z_c * (z[0] ** 2)) + self.Parameters.lam_h_c * 1 / z[2]) * v[0] + (np.sum((self.Problem.station_pow_max * (self.Problem.TOU[0:self.Problem.N_asap] - z[1]))+
                        self.Parameters.lam_z_uc * z[1] ** 2) + self.Parameters.lam_h_uc * 1 / z[2]) * v[1] +  (np.sum(self.Problem.station_pow_max * (self.Problem.TOU[0: self.Problem.N_asap] - 0))) * v[2] + self.Parameters.mu * (conj(v) - np.sum(integ(v)))
        A = np.diag((-1* np.ones((1,3)))[0])
        b = np.zeros((3,1))
        lower_bound= np.zeros((3,1)) 
        beq = np.array([[-(1 - self.Parameters.soft_v_eta)], [1 + self.Parameters.soft_v_eta]])
        model_v = Model("arg_v")
        ## Add all constraints 
        model_v.add_constraint(b[0]-(A[0] @ v).item(0))
        model_v.add_constraint(b[1]-(A[1] @ v).item(0))
        model_v.add_constraint(b[2]-(A[2] @ v).item(0))
        model_v.add_constraint(beq[0] - (Aeq[0]@v).item(0))
        model_v.add_constraint(beq[1] - (Aeq[1]@v).item(0))
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
        x_vars = {(i,1) for i in range(N+N+1)}
        int1 = [x_vars[self.Problem.N_flex + 1 :][i] * (self.Problem.TOU[:self.Problem.N_flex] - z[0])[i] for i in range(N)]
        int2 = [self.Parameters.lam_x * x_vars[self.Problem.N_flex + 1:][i] for i in range(N)]
        # The objective function 
        obj_J = ((cp.sum(int1) + cp.sum(int2) + self.Parameters.lam_z_c * z[0] ** 2) + self.Parameters.lam_h_c * 1 / z[2]) * v[0] + (np.sum((self.Problem.station_pow_max * (self.Problem.TOU[: self.Problem.N_asap] - z[1])) + self.Parameters.lam_z_uc * z[1] ** 2)
               + self.Parameters.lam_h_uc * 1 / z[2]) * v[1] + np.sum(self.Problem.station_pow_max * (self.Problem.TOU[: self.Problem.N_asap ] - 0)) * v[2]

        ## Prepare the data to properly put into the constraints 
        A1L = np.hstack((np.zeros((1, N)), np.array([[-1]])))
        A1R = np.array(np.zeros((1, N)))
        A = np.hstack((A1L , A1R))
        b = -1 * self.Problem.user_SOC_need
        lb1 = np.zeros((N + 1, 1))
        lb1[-1] = self.Problem.user_SOC_need
        lb2 = self.Problem.station_pow_min * np.ones((N,1))
        lb = np.vstack((lb1, lb2))
        ub1 = np.ones((N + 1, 1)) # soc
        ub2 = self.Problem.station_pow_max * np.ones((N,1)) # power
        ub = np.vstack((ub1, ub2))

        # equality constraints - system dynamics
        C1L = np.hstack((np.array([1]).reshape(1,1) ,np.zeros((1, N)))) # initial soc for the first element
        C1R = np.zeros((1, N))

        C2L = np.hstack((np.diag((-1 * np.ones((1, N)))[0]),np.zeros((N, 1)))) + np.hstack((np.zeros((N, 1)), np.diag(np.ones((1, N))[0])))
        C2R = np.diag((-1 * 1.47 * np.ones((1, N)))[0])

        d1 = self.Problem.user_SOC_init
        d2 = np.zeros((N, 1))
        Aeq = np.vstack((np.hstack((C1L, C1R)), np.hstack((C2L, C2R))))
        beq = np.vstack((d1, d2))
        model_x.add_constraint(A @ x_vars <= b)
        model_x.add_constraint(Aeq @ x_vars == beq)
        model_x.add_constraint(lb <= x_vars)
        model_x.add_constraint(x <= ub)
        return model_x.minimize(obj_J)
    def arg_z(self,x,v):
        if sum(v) < 0 | (sum(v) < 1 - self.Parameters.soft_v_eta) | (sum(v) > 1 + self.Parameters.soft_v_eta):
                raise ValueError('[ ERROR] invalid {0}'.format(v))
        ## Create a Model Z that will minimize the objective function 
        model_z = Model("arg_z")

        ## Put the data into proper format for the objective function 
        int0 = lambda z: np.sum([self.Problem.THETA[0][i] * z[i] for i in range(4)])
        int1 = lambda z: np.sum([self.Problem.THETA[1][i] * z[i] for i in range(4)])
        int2 = lambda z: np.sum([self.Problem.THETA[2][i] * z[i] for i in range(4)])
        int3 = lambda z: np.sum([z[i] * (self.Problem.THETA.T @ v)[i] for i in range(4)])
        lse = lambda z: np.log(np.exp(int0(z)) + np.exp(int1(z)) + np.exp(int2(z)))
        J = lambda z: (((np.sum((np.multiply(x[self.Problem.N_flex + 1 :], (self.Problem.TOU[0: self.Problem.N_flex] - z[1]))) +
                         np.multiply(self.Parameters.lam_x, x[self.Problem.N_flex + 1:])) +
                  self.Parameters.lam_z_c * z[0] ** 2) +
                 self.Parameters.lam_h_c * 1 / z[2]) * v[0] + ((np.sum((self.Problem.station_pow_max * (self.Problem.TOU[: self.Problem.N_asap ] - z[1]))+self.Parameters.lam_z_uc * z[1] ** 2) + self.Parameters.lam_h_uc * 1 / z[2]) * v[1])
            + (np.sum(self.Problem.station_pow_max * (self.Problem.TOU[: self.Problem.N_asap])) * v[2]) \
            + self.Parameters.mu * (lse(z) - int3(z)))

        A = [0, 0, 0, 1]
        beq = 1
        constraint1 = {'type': 'ineq', 'fun': lambda z: -(A[1] * z[1] + A[0] * z[0] + A[2] * z[2] + A[3] * z[3]) + b}
        constraint2 = {'type': 'eq', 'fun': lambda z: Aeq[0] * z[0] + Aeq[1] * z[1] + Aeq[2] * z[2] + Aeq[3] * z[3] - beq}
        model_z.add_constraint(contraint1)
        model_z.add_constraint(constraint2)
        model_z.add_constraint(((lb[0], ub[0]), (lb[1], ub[1]), (lb[2], ub[2]), (lb[3], ub[3])))
        return model_z.minimize(J)

