import numpy as np 
import ase
from ase.io.trajectory import Trajectory
from ase.io import read
from ase.io.xyz import write_xyz
from ase import Atoms
from esteem.wrappers.orca import ORCAWrapper

atoms = read("gs_PBE0/xyz/ox.xyz")

owrap_gs = ORCAWrapper()
owrap_es1 = ORCAWrapper()

owrap_gs.setup(nprocs=8,maxcore=4000,orca_cmd="/storage/nanosim/orca5/orca")
owrap_es1.setup(nprocs=8,maxcore=4000,orca_cmd="/storage/nanosim/orca5/orca")

state=["gs","es1"]
calculator=["a","b","c","d","e"]
calc_params_orca_gs = {'basis': 'def2-TZVP', 'target': 0, 'func': 'PBE0', 'disp': True, 'solvent':None}
calc_params_orca_es1 = {'basis': 'def2-TZVP', 'target': 1, 'func': 'PBE0', 'disp': True, 'solvent':None}

angle=[]


for i in range(60, 142, 2):
    #print(i)
    angle.append(i)

    atoms.set_angle(1,0,2,angle=i)
    
    gs = owrap_gs.singlepoint(atoms,'gs',calc_params=calc_params_orca_gs)
    es1 = owrap_es1.singlepoint(atoms,'es1',calc_params=calc_params_orca_es1)

    print("The Ground State energy for angle", i,"is:", gs[0], sep=' ')
    print("The First Excited State energy for angle", i,"is:", es1[0], sep=' ')