import numpy as np 
import ase
from ase.io.trajectory import Trajectory
from ase.io import read
from ase import Atoms
from esteem.wrappers.mace import MACEWrapper

mwrap = MACEWrapper()

atoms = read("gs_PBE0/xyz/ox.xyz")
rand_seed = {'a': 1234,'b': 2345, 'c': 3456, 'd': 4567, 'e': 5678}
suffix = {'rattled_50x6'+rs:rand_seed[rs] for rs in rand_seed}
calc_params_mace = {'target': [0,1],'calc_prefix': "",
                    'calc_suffix': suffix,'calc_seed': "ox", 
                    'calc_dir_suffix': 'rattled'}

state=["gs","es1"]
calculator=["a","b","c","d","e"]
angle =[]
energy = []
force = []
dipole = []

def array_split(energies):
    energy_gs = []
    energy_es1 = []
    for i in range(10):
        if i % 2 == 0:
            energy_gs.append(energies[i])
        else:
            energy_es1.append(energies[i])
    return energy_gs, energy_es1

my_lists = {}

for i in range(60, 122, 2):
    #print(i)
    list_name = f"angle_{i}"
    my_lists[list_name + '_gs'] = []
    my_lists[list_name + '_es1'] = []
    atoms.set_angle(2, 1, 0, angle=i)
    energies, forces, dipoles, calc = mwrap.singlepoint(atoms, 'ox', calc_params=calc_params_mace, forces=True,
                                                        dipole=True)

    # Split energies into ground state and excited state 1
    energy_gs, energy_es1 = array_split(energies)

    # Append the split energies to the dynamically named lists
    my_lists[list_name + '_gs'].extend(energy_gs)
    my_lists[list_name + '_es1'].extend(energy_es1)

# Accessing and printing the dynamically named lists
for key, value in my_lists.items():
    print(f"{key}: {value}")