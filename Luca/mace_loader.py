
from esteem.wrappers.mace import MACEWrapper

mwrap = MACEWrapper()
from ase.io import read
atoms = read("gs_PBE0/xyz/ox.xyz")
rand_seed = {'a': 1234,'b': 2345, 'c': 3456, 'd': 4567, 'e': 5678}
suffix = {'rattled_50x1'+rs:rand_seed[rs] for rs in rand_seed}
calc_params_mace = {'target': [0,1],'calc_prefix': "",
                    'calc_suffix': suffix,'calc_seed': "ox", 
                    'calc_dir_suffix': 'rattled'}

energies, forces, dipoles, calc = mwrap.singlepoint(atoms,'ox',calc_params=calc_params_mace,forces=True,dipole=True)

print(energies)
print(forces)
print(dipoles)
