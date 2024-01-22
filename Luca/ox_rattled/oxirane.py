#ESTEEM script for mycosporine-like amino acid investigations

#!/usr/bin/env python
# coding: utf-8

# 
from esteem import drivers
from esteem import parallel
from copy import deepcopy

# Get default arguments for each script, then edit below
from esteem.tasks.solutes import SolutesTask
from esteem.tasks.solvate import SolvateTask
from esteem.tasks.clusters import ClustersTask
from esteem.tasks.spectra import SpectraTask
from esteem.tasks.ml_training import MLTrainingTask
from esteem.tasks.qmd_trajectories import QMDTrajTask
from esteem.tasks.ml_trajectories import MLTrajTask
from esteem.tasks.ml_testing import MLTestingTask
from esteem.active_learning import create_clusters_tasks,create_mltrain_tasks
from esteem.active_learning import create_mltraj_tasks,create_mltest_tasks
from esteem.active_learning import create_spectra_tasks,add_trajectories

# Set up Wrappers to run underlying packages
from esteem.wrappers import orca, amber, specpycode
from esteem.trajectories import get_trajectory_list

from esteem.wrappers.mace import MACEWrapper
mace_wrapper = MACEWrapper()

solutes_task = SolutesTask()
solvate_task = SolvateTask()
clusters_task = ClustersTask()
spectra_task = SpectraTask()
qmd_task = QMDTrajTask()
train_task = MLTrainingTask()
mltraj_task = MLTrajTask()
test_task = MLTestingTask()

# List of all solvents to use
all_solvents = {None: None}

# List of all solutes to use
all_solutes = {"ox": "oxirane"}
max_iter = 1

# Make a dictionary containing calculators for each Active Learning iteration
train_calcs = []
for tcsuff in ['u']:
    for g in range(max_iter):
        train_calcs.append(f'MACEac{g}{tcsuff}')

# Make a list of random seeds for active learning committee calculator
rand_seed = {'a': 1234,'b': 2345, 'c': 3456, 'd': 4567, 'e': 5678}

# Global Settings
xc_funcs = ['PBE0']
which_func = 0
basis = 'def2-TZVP'
targets = {0: "gs", 1: "es1"}

# Settings for Solutes task
solutes_task.wrapper = orca.ORCAWrapper()
solutes_task.wrapper.setup(nprocs=8,maxcore=4000,orca_cmd="/storage/nanosim/orca5/orca")
solutes_task.script_settings = parallel.get_default_script_settings(solutes_task.wrapper)
solutes_task.nroots = 0
solutes_task.rotate = False
solutes_task.calc_params = None
all_solutes_tasks = {}
for func in xc_funcs:
    for target in targets:
        solutes_task.solvent = None
        solutes_task.basis = basis
        solutes_task.disp = True
        solutes_task.vibrations  = True
        solutes_task.func  = func
        solutes_task.target = target
        solutes_task.directory = f'{targets[target]}_{func}'
        all_solutes_tasks[solutes_task.directory] = deepcopy(solutes_task)

# Settings for Solvate task
solvate_task.wrapper = amber.AmberWrapper(nprocs=8)
solvate_task.script_settings = parallel.get_default_script_settings(solvate_task.wrapper)
solvate_task.script_settings['ntask'] = 8
solvate_task.script_settings['ncpu'] = 1
solvate_task.temp        = 300 # Kelvin
solvate_task.boxsize     = 20  # Angstrom
solvate_task.md_geom_prefix = f"gs_{xc_funcs[which_func]}"
solvate_task.ewaldcut = 12.0 # Angstrom
solvate_task.nheat = 10000
solvate_task.ndens = 100000
solvate_task.nequil = 50000
solvate_task.nsteps = 2000
solvate_task.nsnaps = 2500
solvate_task.use_acpype = True
solvate_task.wrapper.acpype_atom_type = 'gaff'
all_solvate_tasks = {}
MDtemps = [300]
for T in MDtemps:
    for pref in ['','s']: # normal and small boxes
        solvate_task.temp = T
        solvate_task.boxsize = 20 if pref!='s' else 12 #13
        solvate_task.ewaldcut = 12.0 if pref!='s' else 9 #9
        solvate_task.md_suffix = f'{pref}md{T}'
        all_solvate_tasks[solvate_task.md_suffix] = deepcopy(solvate_task)

# Generic settings for Clusters task
truth = 'orca'
meth = ''
clusters_task.output   = truth
clusters_task.func     = xc_funcs[which_func]
clusters_task.basis    = basis
clusters_task.disp     = True
clusters_task.boxsize  = None
clusters_task.nroots   = 0
clusters_task.calc_forces = True
clusters_task.radius   = None
all_clusters_tasks = {}

dft_wrapper_8 = orca.ORCAWrapper()
dft_wrapper_8.setup(nprocs=8,maxcore=3000,orca_cmd="/storage/nanosim/orca5/orca")
dft_wrapper_24 = orca.ORCAWrapper()
dft_wrapper_24.setup(nprocs=24,maxcore=3000,orca_cmd="/storage/nanosim/orca5/orca")
dft_wrapper_24.set_careful_scf(sthresh=3e-8)
dft_script_settings_8 = parallel.get_default_script_settings(dft_wrapper_8)
dft_script_settings_8['ntask'] = 1
dft_script_settings_8['ncpu'] = 8
dft_script_settings_8['partition'] = 'nanosimd'
dft_script_settings_8['logdir'] = 'orca'
dft_script_settings_24 = parallel.get_default_script_settings(dft_wrapper_24)
dft_script_settings_24['ntask'] = 1
dft_script_settings_24['ncpu'] = 24
dft_script_settings_24['partition'] = 'nanosim'
dft_script_settings_24['logdir'] = 'orca'
dft_script_settings_24['declarations'] = dft_script_settings_24['declarations'].replace("SMAX=1","SMAX=5")
ml_wrapper = mace_wrapper
ml_wrapper.script_settings = parallel.get_default_script_settings(ml_wrapper)

# Clusters tasks to run Distorted inputs
clusters_task.wrapper = dft_wrapper_8
clusters_task.script_settings = dft_script_settings_8
clusters_task.md_prefix = 'distorted'
clusters_task.md_suffix = 'distort_gs_A_nocalc'
clusters_task.exc_suffix = 'distort'
clusters_task.ref_mol_dir = None
all_clusters_tasks[clusters_task.exc_suffix] = deepcopy(clusters_task)

# Clusters tasks to run Rattled inputs
clusters_task.md_prefix = 'rattled'
clusters_task.exc_suffix = 'rattled'
clusters_task.ref_mol_dir = None
clusters_task.target = [0,1]
for w in get_trajectory_list(7): 
    clusters_task.md_suffix = f'rattled_gs_{w}_nocalc'
    clusters_task.which_traj = w
    all_clusters_tasks[f'{clusters_task.exc_suffix}_{w}'] = deepcopy(clusters_task)

# Clusters tasks to run output of solvate MD
MDchar = ['R']
cluster_radii = {'A': 0.0, 'B': 2.5, 'C': 5.0, 'D': 7.5, 'E': 10.0, 'U': None}
solv_rad = {}
solv_rad[cluster_radii['B']] = {solv: 2.5 for solv in all_solvents}
for i,MDtemp in enumerate(MDtemps):
    for cr in cluster_radii:
        r = cluster_radii[cr]
        clusters_task.md_suffix = 'solv'
        clusters_task.wrapper = dft_wrapper_8 if r==0.0 else dft_wrapper_24
        clusters_task.script_settings = dft_script_settings_8 if r==0.0 else dft_script_settings_24
        MDstr = MDchar[i]
        clusters_task.md_prefix = f'{{solu}}_{{solv}}_smd{MDtemp}'
        clusters_task.exc_suffix = f'solv{MDstr}{r}'
        if r in solv_rad:
            clusters_task.radius = solv_rad[r]
        else:
            clusters_task.radius = r
        for w in ['A','B']:
            clusters_task.target = 0
            clusters_task.nroots = 0
            clusters_task.output = 'orca'
            clusters_task.which_traj = w
            clusters_task.max_snapshots = 200 if w == 'A' else 2040
            clusters_task.min_snapshots = 0 if w == 'A' else 2000
            all_clusters_tasks[clusters_task.exc_suffix+'_'+w] = deepcopy(clusters_task)

# Main active learning loop clusters_tasks
clusters_task.max_snapshots = 100
clusters_task.min_snapshots = 0
clusters_task.valid_snapshots = 20
clusters_task.radius = None
clusters_task.repeat_without_solute = False
clusters_task.wrapper = dft_wrapper_8
clusters_task.subset_selection_nmax = 100
clusters_task.subset_selection_min_spacing = 20
clusters_task.subset_selection_bias_beta = 10000
seed = '{solu}'
traj_suffix = 'mlclus'
md_suffix = 'mldyn_recalc'
md_dir_suffix = 'mldyn'
all_clusters_tasks.update(create_clusters_tasks(clusters_task,train_calcs,seed,traj_suffix,
                                                md_suffix,md_dir_suffix,targets,rand_seed,meth,truth))

# Quantum Molecular Dynamics
from ase.units import fs,AUT
all_qmd_tasks = {}
qmd_task.wrapper = orca.ORCAWrapper()
qmd_task.script_settings = parallel.get_default_script_settings(qmd_task.wrapper)
qmd_task.wrapper.setup(nprocs=4,maxcore=3200,orca_cmd="/storage/nanosim/orca5/orca") # WAS 32
qmd_task.temp         = 800
qmd_task.nsnap        = 20
qmd_task.nequil       = 10
qmd_task.qmd_steps    = 20
qmd_task.qmd_timestep = 0.5*fs
qmd_task.ntraj        = 32
qmd_task.constraints  = None
qmd_task.func         = xc_funcs[which_func]
qmd_task.basis        = basis
qmd_task.disp         = True
qmd_task.solvent      = None #'{solu}'

# add ground and excited state tasks
for target in targets:
    qmd_task.target       = [target,1-target]
    qmd_task.traj_suffix  = f'{targets[target]}_qmd{qmd_task.temp}'
    qmd_task.geom_prefix  = f'{targets[target]}_{qmd_task.func}'
    all_qmd_tasks[qmd_task.traj_suffix] = deepcopy(qmd_task)

# add a version for use in the atoms task (no excitations)
qmd_task.ref_mol_dir  = f'gs_{xc_funcs[which_func]}'
qmd_task.target = 0
qmd_task.traj_suffix = 'orca'
qmd_task.solvent = '{solv}'
all_qmd_tasks[qmd_task.traj_suffix] = deepcopy(qmd_task)

# Generic Training task
train_task.wrapper = mace_wrapper
train_task.script_settings = parallel.get_default_script_settings(train_task.wrapper)
train_task.cutoff = 5.0
train_task.restart = False
train_task.ntraj = 270
train_task.reset_loss = False
# Set MACE-specific parameters
train_task.wrapper.train_args['batch_size'] = 16
train_task.wrapper.train_args['valid_batch_size'] = 16
train_task.wrapper.train_args['dipole_key'] = 'dipole'
train_task.wrapper.train_args['loss'] = 'energy_forces_dipole'
train_task.wrapper.train_args['max_ell'] = 3
train_task.wrapper.train_args['num_interactions'] = 2
train_task.wrapper.train_args['max_num_epochs'] = 1000
train_task.wrapper.train_args['swa'] = True
train_task.wrapper.train_args['start_swa'] = 900
train_task.wrapper.train_args['model'] = 'EnergyDipolesMACE'
train_task.wrapper.train_args['hidden_irreps'] = '32x0e + 32x1o + 32x2e'
train_task.wrapper.train_args['error_table'] = 'EnergyDipoleRMSE'
train_task.wrapper.train_args['restart_latest'] = True

# Create active learning mltrain tasks
traj_suffixes = []
dir_suffixes = {}
ntraj = {}
# Make a list of the input trajectories (just "rattled" here)
if True:
    traj_suffixes.append("rattled")
    dir_suffixes["rattled"] = "rattled"
    ntraj[targets[0],"rattled"] = 1
    ntraj[targets[1],"rattled"] = 0
else:
    traj_suffixes.append("solvR2.5")
    dir_suffixes["solvR2.5"] = "solvR2.5"
    ntraj[targets[0],"solvR2.5"] = 1
    ntraj[targets[1],"solvR2.5"] = 0

# Option to add distorted trajectories
if False:
    traj_suffixes.append("distorted")
    dir_suffixes["distorted"] = "distorted"
    ntraj[targets[0],"distorted"] = 1
    ntraj[targets[1],"distorted"] = 0

# Directory suffix for iterative learning clusters
iter_dir_suffixes = ["mlclus"]
seeds=["{solu}"]
all_mltrain_tasks = {}
all_mltrain_tasks.update(create_mltrain_tasks(train_task,train_calcs,seeds,targets,rand_seed,meth,truth,traj_suffixes,dir_suffixes,ntraj,iter_dir_suffixes,delta_epochs=200,separate_valid=True))

# Non-active learning calculators (rattled configs only)
train_task.wrapper.train_args['max_num_epochs'] = 1000
train_task.wrapper.train_args['swa'] = True
train_task.wrapper.train_args['start_swa'] = 800
rattled_calcs = [f'rattled_50x{i}' for i in range(1,11)]
seeds = ["{solu}"]
for target in targets:
    targstr = targets[target]
    for t in rattled_calcs:
        train_task.reset_loss = False
        targstr = targets[target]
        train_task.traj_suffix = truth
        train_task.calc_dir_suffix = 'rattled'
        train_task.target = target
        train_task.which_trajs = []
        # Set up links to trajectories
        train_task.traj_links = {}
        # Add n trajectories
        ntraj[targstr,"rattled"] = int(t[11:])
        
        ntraj[targstr,"rattled"] = 0

        for rs in rand_seed:
            train_task.wrapper.train_args['seed'] = rand_seed[rs] # MACE specific
            train_task.calc_suffix = f'{t}{rs}'
            all_mltrain_tasks[f'{targstr}_{t}{rs}'] = deepcopy(train_task)

# Settings for ML Molecular Dynamics
md_wrapper = mace_wrapper
snap_wrapper = deepcopy(mace_wrapper)
mltraj_task.script_settings = parallel.get_default_script_settings(mltraj_task.wrapper)
mltraj_task.temp        = 300 # Kelvin
mltraj_task.thermostat = "NPT"
mltraj_task.md_timestep = {'MD': 0.5*fs, 'EQ': 0.5*fs}
mltraj_task.md_friction = {'MD': 0.002, 'EQ': 0.05}  # For Langevin dynamics
#mltraj_task.md_friction = {'MD': 15*fs, 'EQ': 3*fs} # For NPT dynamics
# Number of MD steps per snapshot
mltraj_task.md_steps    = 10
# Number of different trajectories (ABCDE...) and list of trajectories -changed from 5 to 7
mltraj_task.ntraj       = 7
mltraj_task.which_trajs = get_trajectory_list(mltraj_task.ntraj)
# Number of snapshots per trajectory (each one md_steps in length)
# Here we do 10000 steps, saving every 10th one
mltraj_task.nsnap       = 1000
# Number of snapshots of equilibration (each one md_steps in length)
mltraj_task.nequil      = 100
mltraj_task.constraints = None
mltraj_task.store_full_traj = False
# Once the trajectory is finished, remove solvent molecules beyond 2.5A of solute
mltraj_task.carve_trajectory_radius = 2.5
# Recalculate with the snap_wrapper after carving, to get error estimate
mltraj_task.recalculate_carved_traj = True
mltraj_task.continuation = True
mltraj_task.dynamics = None
mltraj_task.geom_prefix  = f'gs_{xc_funcs[which_func]}'
mltraj_task.md_init_traj_link = f"{{solu}}_{{solv}}_smd300/{{solu}}_{{solv}}_solv.traj"
mltraj_task.calc_seed = f"{list(all_solutes)[0]}_{{solv}}"

# Create active learning mltraj tasks
all_mltraj_tasks = {}
all_mltraj_tasks.update(create_mltraj_tasks(mltraj_task,train_calcs,targets,rand_seed,meth,md_wrapper,snap_wrapper))
mltraj_task.nsnap       = 80000
mltraj_task.traj_suffix = "specdyn"
all_mltraj_tasks.update(create_mltraj_tasks(mltraj_task,train_calcs,targets,rand_seed,meth,md_wrapper,snap_wrapper))

# Settings for ML_Testing task
test_task.wrapper = mace_wrapper
test_task.script_settings = parallel.get_default_script_settings(test_task.wrapper)
test_task.ntraj = 300
test_task.ref_mol_dir = f'{{targ}}_{xc_funcs[which_func]}'
test_task.calc_seed = f"{list(all_solutes)[0]}_{{solv}}"
all_test_tasks = {}
# Add the rattled_B-J trajectories as extra test sets
ntraj[targets[0],"rattled"] = 1
ntraj[targets[1],"rattled"] = 0
all_train_calcs = rattled_calcs #train_calcs + rattled_calcs
all_test_tasks.update(create_mltest_tasks(test_task,all_train_calcs,seeds,targets,rand_seed,truth,meth,
                                          traj_suffixes,dir_suffixes,iter_dir_suffixes,ntraj))

#test task that tests only on unseen, final calc generation produced trajectories
meth = ''
for target in targets:
        for t in train_calcs:
            for rs in rand_seed:
                test_task.wrapper.train_args['seed'] = rand_seed[rs]
                test_task.calc_suffix = f'{meth}{t}{rs}'
                test_task.plotfile = f'{{solu}}_{test_task.calc_suffix}.png'
                # This test uses the calculator directory from the MLTrain task as the traj location
                test_task.traj_suffix = truth
                test_task.calc_prefix = ""
                test_task.calc_dir_suffix = f'{meth}{t[0:6]}'
                test_task.which_trajs = list('A')
                test_task.traj_prefix = f"{test_task.calc_seed}_{targets[target]}_{meth}{t[0:6]}_test/"
                test_task.target = target
                targstr = targets[target]
                test_task.traj_links = {}
                test_task.which_trajs = []
                if True:
                    p = ''
                    targ0 = "gs"
                    v = 'B'
                    #meth = 'MACE'
                    method = "orca"
                    gen_test = 2
                    sel_char = 'r'
                    test_task.which_trajs = ['B']
                    test_task.traj_prefix = f'{{solu}}_{{solv}}_{targets[target]}_MACEac_test/'
                    for w in test_task.which_trajs:
                        #if w in test_task.calc_suffix:
                         #   continue
                        test_task.traj_links[v+w] = f'{{solu}}_{{solv}}_{targ0}_MACEac_mlclus/{{solu}}_{{solv}}_{targets[target]}_R_{method}_ac{gen_test}{sel_char}.traj'
                    test_task.which_trajs = [v+'B']
                    all_test_tasks[f"{targets[target]}_{meth}{t}{rs}_mltraj_MACEac{gen_test}{sel_char}"] = deepcopy(test_task)
            
# Settings for Spectra script

# Direct spectra from ORCA calculations
all_spectra_tasks = {}
spectra_task.broad = 0.05 # eV
spectra_task.inputformat   = 'orca'
spectra_task.wavelength    = (300,700,1) # nm
spectra_task.warp_origin_prefix = f'{xc_funcs[which_func]}/is_tddft_{{solv}}'
spectra_task.warp_dest_prefix   = f'{xc_funcs[which_func]}/is_tddft_{{solv}}'
spectra_task.warp_inputformat   = 'orca'
spectra_task.warp_files         = '{solu}.out'
spectra_task.warp_scheme        = 'beta'
spectra_task.exc_suffix         = 'IS'
all_spectra_tasks['IS'] = deepcopy(spectra_task)

# SpecPyCode spectra based on MD trajectories from mltraj calculations
spectra_task.inputformat = 'traj'
spectra_task.warp_scheme = None
spectra_task.warp_files  = None
spectra_task.files       = None
spectra_task.wrapper = specpycode.SpecPyCodeWrapper(rootname="spec_test",input_filename="spec_test_input")
spectra_task.start_frame = 1000
spectra_task.max_frames = 5000
spectra_task.ncores = 6
if spectra_task.wrapper is not None:
    spectra_task.wrapper.md_step = mltraj_task.md_timestep['MD']/fs*mltraj_task.md_steps
    spectra_task.wrapper.num_steps = spectra_task.max_frames
    spectra_task.wrapper.temperature = mltraj_task.temp
    spectra_task.wrapper.max_t = spectra_task.max_frames * spectra_task.wrapper.md_step
    spectra_task.wrapper.decay_length = 500
    spectra_task.wrapper.corr_length_3rd = 1000
    spectra_task.wrapper.spectral_window = 1.8 # eV
    spectra_task.wrapper.solvent_model = 'NONE'
    spectra_task.wrapper.num_trajs = 1

#all_spectra_tasks.update(create_spectra_tasks(spectra_task,train_calcs,targets,rand_seed,meth))

import sys
if "ipykernel_launcher.py" in sys.argv[0]:
    raise Exception("Skipping main routine execution as notebook environment detected")

drivers.main(all_solutes, all_solvents,
             all_solutes_tasks=all_solutes_tasks,
             all_solvate_tasks=all_solvate_tasks,
             all_clusters_tasks=all_clusters_tasks,
             all_mltrain_tasks=all_mltrain_tasks,
             all_mltraj_tasks=all_mltraj_tasks,
             all_mltest_tasks=all_test_tasks,
             all_qmd_tasks=all_qmd_tasks,
             all_spectra_tasks=all_spectra_tasks,
             make_script=parallel.make_sbatch)

exit()

