
from ase.io import read

atoms = read("gs_PBE0/xyz/ox.xyz")

from esteem.wrappers.orca import ORCAWrapper

owrap_gs = ORCAWrapper()
owrap_es1 = ORCAWrapper()

owrap_gs.setup(nprocs=8,maxcore=4000,orca_cmd="/storage/nanosim/orca5/orca")
owrap_es1.setup(nprocs=8,maxcore=4000,orca_cmd="/storage/nanosim/orca5/orca")

calc_params_orca_gs = {'basis': 'def2-TZVP', 'target': 0, 'func': 'PBE0', 'disp': True, 'solvent':None}
calc_params_orca_es1 = {'basis': 'def2-TZVP', 'target': 1, 'func': 'PBE0', 'disp': True, 'solvent':None}

a=owrap_gs.singlepoint(atoms,'groundstate',calc_params=calc_params_orca_gs)
b=owrap_es1.singlepoint(atoms,'excitedstate1',calc_params=calc_params_orca_es1)

print(a[0])
print(b[0])