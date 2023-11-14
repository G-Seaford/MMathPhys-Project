from ase.io import read,Trajectory
import numpy as np

opt = read("../es1_PBE0/opt/ox_es1.xyz")
seeds={"A":42,"B":1234,"C":2345,"D":3456,"E":4567,
       "F":56,"G":5678,"H":6789,"I":7890,"J":8901}

for w in ["A","B","C","D","E"]:
    outtrajname = f"ox_rattled_es1_{w}_nocalc.traj"
    traj = Trajectory(outtrajname,"w")
    rng = np.random.RandomState(seeds[w])
    for i in range(0,50):
       rat = opt.copy()
       rat.rattle(0.05,rng=rng)
       print(rat.positions[0])
       traj.write(rat)
    traj.close()

    outfilename=f"ox_rattled_es1_{w}_nocalc.xyz"
    from ase.io.extxyz import write_xyz
    f=open(outfilename,"w")
    write_xyz(f,Trajectory(outtrajname))
    f.close()

