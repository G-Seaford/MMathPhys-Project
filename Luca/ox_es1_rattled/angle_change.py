import ase
from ase.io.trajectory import Trajectory
from ase import Atoms
import numpy as np 


name="ox_es1_orca_merged_rattled_50x1a.traj"

traj = Trajectory(name)
atoms = traj[-1]

#pe=atoms.get_potential_energy()
#print(pe)
#ke=atoms.get_kinetic_energy()
#print(ke)

atomicn=atoms.get_positions()
oxygen=atoms.get_positions()[0]
carbon1=atoms.get_positions()[1]
carbon2=atoms.get_positions()[2]

atomicn2=atoms.get_atomic_numbers()
print(atomicn)
print(atomicn2)
angles=atoms.get_angle(1,0,2)
print(angles)

atoms.set_angle(1,0,2,angle=108)
print(atoms.get_angle(1,0,2))

#atomicn3=atoms.get_potential_energy()

#print(atomicn3)

atoms.get_potential_energy()
#print(angles)

##pes=atoms.get_potential_energies()
##print(pes)