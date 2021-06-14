import os
import numpy as np
import matplotlib.pyplot as plt

folder = "./run6/"
dir_list_file = "./run6/pfolder.txt"
file = "/trajStates.dat"

seeds = np.loadtxt(dir_list_file,dtype=str)
iflag = 0
for s in seeds:
	fpath = folder+s+file

	tim     = np.loadtxt(fpath, usecols=(0,))
	count   = np.loadtxt(fpath, usecols=(4,5))
	if(iflag==0):
		ctot    = count
		iflag 	= 1
	else: 
		ctot    = np.add(ctot,count)
	print(fpath, ctot.shape)

ptot = np.sum(ctot,1) 
print(ptot)
pop1  = (ctot[:,0])/ptot
pop12 = (ctot[:,0]-ctot[:,1])/ptot

eng2au = 4.55634e-6
au2fs = 2.4e-2 
timplot = tim*au2fs
yval = pop1
print(len(timplot),len(yval))

np.savetxt('trajdata_all.dat', np.c_[timplot, yval])
plt.plot(timplot, yval)
plt.xlabel('time (fs)')
plt.ylabel('population1')
# plt.ylabel('population diff')
# plt.ylim([0.0,1.1])
plt.grid()
# plt.show()

plt.savefig("pop.png")
