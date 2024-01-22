###########################################################################################
# The F_IJ Calculator for MACE to perform a BFGS optimisation.
# Authors: Gianluca Seaford, Ovidiu Tirnovan
###########################################################################################

import numpy as np

from ase import io
from ase import Atoms
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.calculators.orca import ORCA
from ase.constraints import FixBondLengths, FixLinearTriatomic
from ase.vibrations import Vibrations
from numpy import cross, eye, dot
from scipy.linalg import expm, norm
from esteem.wrappers.mace import MACEWrapper
from esteem.wrappers.orca import ORCAWrapper
from mace.calculators import MACECalculator






class F_IJ_MACECalculator(MACECalculator):
    
    def __init__(self, alpha, delta, sigma, tolstep, tolgrad):
        '''Initialises the parent class'''
        super().__init__()
        
        '''COnstants for this calculator'''
        self.alpha = 0.5442281591 # eV
        self.delta = 0.02721140795 # eV
        self.sigma = 8 # Dimensionless Constant
        self.tolstep = 2.721140795e-5 # eV
        self.tolgrad = 0.257110524 # eV/Å
        
    @staticmethod
    def array_split(input, len_array):
        '''Function to split arrays into S0 and S1 components'''
        gs = []
        es1 = []
        for i in range(len_array):
            if i % 2 == 0:
                gs.append(input[i])
            else:
                es1.append(input[i])
        return gs, es1
    
    def find_mean(self, atoms, positions):
        '''Calculates the energy array, splits this into S_0 and S_1 contributions and returns these arrays'''
        atoms = Atoms(symbols=atoms.get_chemical_symbols(), positions=positions.reshape((-1, 3)))
        all_energy, _, _, _ = self.mwrap.singlepoint(atoms, 'ox', calc_params=self.calc_params_mace, forces=True, dipole=True)
        energy_gs, energy_es1 = array_split(all_energy)
        energy_gs_mean = np.mean(energy_gs)
        energy_es1_mean = np.mean(energy_es1)
        return energy_gs_mean,energy_es1_mean
    
    def E_IJ(self, energy_gs_mean, energy_es1_mean):
        '''Calculates mean midpoint between the two states '''
        return np.mean(energy_gs_mean + energy_es1_mean)
    
    def Delta_E_IJ(self, energy_gs_mean, energy_es1_mean):
        '''Calculates mean energy difference between the two states '''
        return energy_gs_mean - energy_es1_mean
    
    def G_IJ(self, delta_E_IJ, alpha):
        '''The Penalty function used in the calculation of F'''
        delta_E_IJ_squared = np.square(delta_E_IJ)
        return np.divide(delta_E_IJ_squared, (np.array(delta_E_IJ) + alpha))
        
    def F_IJ(self, sigma, alpha, atoms, positions):
        '''Takes input parameters, calls relevant functions and outputs F_IJ'''
        energy_gs, energy_es1 = self.find_mean(atoms, positions)
        delta_IJ = self.Delta_E_IJ(energy_gs, energy_es1)
        f_IJ = self.E_IJ(energy_gs, energy_es1) + self.sigma * self.G_IJ(delta_IJ, alpha)
        return f_IJ
    
    def dF_IJ(self, alpha, delta_E_IJ, f_I, f_J):
        '''Calculates the derivatives of F_IJ'''
        delta_E_IJ_squared = np.square(delta_E_IJ)
        dF_IJ = 1/2 *(f_i + f_j) + self.sigma * (delta_E_IJ_squared + 2*self.alpha*delta_E_IJ)/np.square(delta_E_IJ + self.alpha) * (f_I - f_J)
        return dF_IJ
        
    def calculate(self, atoms, init_config = None):
        if init_config is None:
            init_config = atoms.get_positions()
        MACECalculator.calculate()
        
        
        
        #Recombine