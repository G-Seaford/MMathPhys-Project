import numpy as np 
import ase
from ase.io.trajectory import Trajectory
from ase.io import read
from ase import Atoms
from esteem.wrappers.mace import MACEWrapper

mwrap = MACEWrapper()

# Defines file names for ease of use #
fname_gs_a="ox_gs_A_orca.traj"
fname_es1_a="ox_es1_A_orca.traj"

# Gets trajectories from files and assigns them to variables #
t_0 = Trajectory(fname_gs_a)
t_1 = Trajectory(fname_es1_a)

# Gets the positions of the atoms from the trajectory variables #
atoms_gs = t_0[-1]
atoms_es1 = t_1[-1]

# Gets potential energy and atom locations of gs and es1 configurations #
pe_gs = atoms_gs.get_potential_energy()
pe_es1 = atoms_es1.get_potential_energy()
locations_gs = atoms_gs.get_positions()
locations_es1 = atoms_es1.get_positions()

# Finds locations of the Oxygen and Carbon atoms in the ground state #
oxygen_gs = locations_gs[0]
carbon_1_gs = locations_gs[1]
carbon_2_gs = locations_gs[2]
coc_angle_gs = atoms_gs.get_angle(1,0,2)
print("Ground State Atom Configurations")
print("Oxygen:",oxygen_gs, sep=" ") 
print("Carbon 1:",carbon_1_gs, sep=" ") 
print("Carbon 2:",carbon_2_gs, sep=" ")
print("C-O-C Angle:",coc_angle_gs, sep=" ")
print(" ")

# Finds locations of the Oxygen and Carbon atoms in the first excited state #
oxygen_es1 = locations_es1[0]
carbon_1_es1 = locations_es1[1]
carbon_2_es1 = locations_es1[2]
coc_angle_es1 = atoms_es1.get_angle(1,0,2)
print("First Excited State Atom Configurations")
print("Oxygen:",oxygen_es1, sep=" ") 
print("Carbon 1:",carbon_1_es1, sep=" ") 
print("Carbon 2:",carbon_2_es1, sep=" ")
print("C-O-C Angle:",coc_angle_es1, sep=" ")
print(" ")


#Loops over angles for configurations #
gs_config = []
es1_config = []
for i in range(60, 145, 5):
    # Perturbs the Ground State positions #
    gs = atoms_gs.set_angle(1,0,2,angle=i)
    gs=atoms_gs.get_positions()
    
    # Perturbs the First Excited State Parameters #
    es1 = atoms_es1.set_angle(1,0,2,angle=i)
    es1 = atoms_es1.get_positions()
    
    # Outputs the perturbed arrays into a list for later use #
    gs_config.append(gs) 
    es1_config.append(es1)

atoms = read("gs_PBE0/xyz/ox.xyz")
rand_seed = {'a': 1234,'b': 2345, 'c': 3456, 'd': 4567, 'e': 5678}
suffix = {'rattled_50x1'+rs:rand_seed[rs] for rs in rand_seed}
calc_params_mace = {'target': [0,1],'calc_prefix': "",
                    'calc_suffix': suffix,'calc_seed': "ox", 
                    'calc_dir_suffix': 'rattled'}

energies, forces, dipoles, calc = mwrap.singlepoint(atoms,'ox',calc_params=calc_params_mace,forces=True,dipole=True)

print(energies)
print(forces)
