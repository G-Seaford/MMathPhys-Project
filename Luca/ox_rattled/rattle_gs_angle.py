from ase.io import read,Trajectory
import numpy as np

opt = read("../gs_PBE0/opt/ox.xyz")
seeds={"A":42,"B":1234,"C":2345,"D":3456,"E":4567,
       "F":56,"G":5678,"H":6789,"I":7890,"J":8901}
trajs = ["A","B","C","D","E","F","G"]

for j,theta in enumerate(range(60, 130, 10)):
    opt.set_angle(2,1,0,angle=theta)
    w = trajs[j]
    outtrajname = f"ox_rattled_gs_{w}_nocalc.traj"
    traj = Trajectory(outtrajname,"w") # Check with Nick if w needed instead of a #
    rng = np.random.RandomState(seeds[w])
    nconfigs = 50 if j==0 else 10
    for i in range(0,nconfigs):
        rat = opt.copy()
        rat.rattle(0.05,rng=rng)
        print(rat.positions[0])
        traj.write(rat)
    traj.close()

    outfilename=f"ox_rattled_gs_{w}_nocalc.xyz"
    from ase.io.extxyz import write_xyz
    f=open(outfilename,"w")
    write_xyz(f,Trajectory(outtrajname))
    f.close()
