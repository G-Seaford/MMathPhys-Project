import scipy as sp
import numpy as np
from ase import io
from ase import Atoms
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.vibrations import Vibrations
from ase.calculators.orca import ORCA
from ase.constraints import FixBondLengths, FixLinearTriatomic
from numpy import cross, eye, dot
from scipy.linalg import expm, norm
from esteem.wrappers.mace import MACEWrapper
from esteem.wrappers.orca import ORCAWrapper


mwrap = MACEWrapper()
atoms = read("gs_PBE0/xyz/ox.xyz")
rand_seed = {'a': 1234,'b': 2345, 'c': 3456, 'd': 4567, 'e': 5678}
suffix = {'rattled_50x5'+rs:rand_seed[rs] for rs in rand_seed}
calc_params_mace = {'target': [0,1],'calc_prefix': "",
                    'calc_suffix': suffix,'calc_seed': "ox", 
                    'calc_dir_suffix': 'rattled'}
def array_split(energies):
    energy_gs = []
    energy_es1 = []
    for i in range(10):
        if i % 2 == 0:
            energy_gs.append(energies[i])
        else:
            energy_es1.append(energies[i])
    return energy_gs, energy_es1

atoms_pos = atoms.get_positions()

calc_params_mace = {'target': [0,1],'calc_prefix': "",
                    'calc_suffix': suffix,'calc_seed': "ox", 
                    'calc_dir_suffix': 'rattled'}

##initial values (make sure units are consistent with ORCA)

##############################################################################################################################

alpha = 0.5442281591 # eV
delta = 0.02721140795 # eV
sigma = 11 # Dimensionless Constant
tolstep = 2.721140795e-5 # eV
tolgrad = 0.257110524 # eV/Å

def EIJ(atoms, positions):
    "Function that gives E_IJ, i.e the average of the 2 differences, a mean is taken for all 5 calculators"
    atoms = Atoms(symbols=atoms.get_chemical_symbols(), positions=positions.reshape((-1, 3)))
    energies, _, _, _ = mwrap.singlepoint(atoms, 'ox', calc_params=calc_params_mace, forces=True, dipole=True)
    energy_gs, energy_es1 = array_split(energies)
    energygs=np.mean(energy_gs)
    energyes1=np.mean(energy_es1)
    print("GS:",energygs,"ES1:",energyes1)
    return np.mean(energygs + energyes1)  # Calculate the mean before returning


def deltaEIJ(atoms, positions):
    "function that gives the difference between the 2 energy values, mean taken for all 5 calculators"
    atoms = Atoms(symbols=atoms.get_chemical_symbols(), positions=positions.reshape((-1, 3)))
    atoms.set_calculator(mwrap)
    energies, _, _, _ = mwrap.singlepoint(atoms, 'ox', calc_params=calc_params_mace, forces=True, dipole=True)
    energy_gs, energy_es1 = array_split(energies)

    energygs = np.mean(energy_gs)
    energyes1 = np.mean(energy_es1)

    # Compute the difference between the means
    delta_energy = energyes1 - energygs

    return delta_energy

def GIJ(deltaEIJ, alpha):
    "G_IJ function"
    deltaEIJ_squared = np.square(deltaEIJ)
    return np.divide(deltaEIJ_squared, (np.array(deltaEIJ) + alpha))



def FIJ(positions, sigma, alpha):
    "Final F_IJ function as requested, function which is iterated by algorithm"
    deltasEIJ = deltaEIJ(atoms, positions)
    result = EIJ(atoms, positions) + sigma * GIJ(deltasEIJ, alpha)
    #positions = positions.reshape((-1, 3))  # Reshape to 2D array

    #print("Coordinates:", positions)

    return result



def FORCE(atoms):
    
    #fix code in terms of atoms at the particular point
    
    energies, forces, dipoles, calc = mwrap.singlepoint(atoms, 'ox', calc_params=calc_params_mace, forces=True,
                                                        dipole=True)
    forces_gs, forces_es1 = array_split(energies)
    return forces_gs,forces_es1


print(atoms_pos)

#atoms_pos = atoms.get_positions().flatten()

result = sp.optimize.minimize(FIJ, atoms_pos, args=(sigma, alpha),method='L-BFGS-B')
