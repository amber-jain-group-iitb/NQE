# Time-Averaging over trajectories run on different cores
import numpy as np
import matplotlib.pyplot as plt 

folder          = "./run1/"
folders = np.loadtxt(folder+"pfolder.txt", dtype='str')
p       = np.loadtxt("input0.txt")

ttraj = int(p[1])*int(p[2])                 # ncore*ntraj
ntsteps = np.ceil(p[4]/p[3]).astype('int')  # ttime/dtc
ns = np.ceil(ntsteps/p[8]).astype('int')    # ntsteps/nwrt_traj
print(ttraj,ntsteps, ns)

lpop =10
popsums = np.zeros(ns)
for f in folders:
	fpath = folder + f + "/fort.100"
	p = np.loadtxt(fpath, dtype='float')
	tsteps, psum = zip(*p)
	popsums = np.add(popsums, psum)
popsums_avg = popsums/ttraj
longpop = np.mean( popsums_avg[-lpop:] )

np.savetxt('pops.dat', np.vstack((tsteps, popsums_avg)).T, fmt='%18.10e')
f = open("longpops.txt", "w")
f.write(str(longpop))
f.close()

