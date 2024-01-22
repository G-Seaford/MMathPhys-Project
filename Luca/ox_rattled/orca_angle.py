import numpy as np 
import ase
from ase.io.trajectory import Trajectory
from ase.io import read
from ase import Atoms
from esteem.wrappers.orca import ORCAWrapper

owrap = ORCAWrapper()
owrap.nprocs = 8

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
coc_angle_gs = atoms_gs.get_angle(2,1,0)
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
coc_angle_es1 = atoms_es1.get_angle(2,1,0)
print("First Excited State Atom Configurations")
print("Oxygen:",oxygen_es1, sep=" ") 
print("Carbon 1:",carbon_1_es1, sep=" ") 
print("Carbon 2:",carbon_2_es1, sep=" ")
print("C-O-C Angle:",coc_angle_es1, sep=" ")
print(" ")


#Loops over angles for configurations #
gs_config = []
es1_config = []
for i in range(60, 142, 2):
    # Perturbs the Ground State positions #
    gs = atoms_gs.set_angle(2,1,0,angle=i)
    gs=atoms_gs.get_positions()
    
    # Perturbs the First Excited State Parameters #
    es1 = atoms_es1.set_angle(2,1,0,angle=i)
    es1 = atoms_es1.get_positions()
    
    # Outputs the perturbed arrays into a list for later use #
    gs_config.append(gs) 
    es1_config.append(es1)
    print(i, gs_config[i], es1_config[i], sep=",")

#atoms= gs_config[1]

#print("Ground State:",gs_config,sep=" ")
#print("First Excited State:",es1_config,sep=" ")


#calc_params_orca = {'basis': 'def2-TZVP', 'target': 0, 'func': 'PBE0', 'disp': True, 'solvent':None}

#owrap.singlepoint(atoms,'test',calc_params=calc_params_orca)


