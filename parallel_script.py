import subprocess
import numpy as np
import os.path

compiled_fname 	= "a.out"
job_script     	= "sub.sh"
input_0 		= "input0.txt"
input_i			= "input.txt"
folder			= "run1/"
dir_list		= "pfolder.txt"

params	= ["inparam_tk.txt", "inparam_vc.txt","inparam_l1.txt","inparam_et.txt"]

p = np.loadtxt(input_0)
ncore = int(p[1])
#Input text files
for seed in range(1,ncore+1):
	np.savetxt(input_i, p)
	file = open(input_i, 'a')
	file.write(str(seed) +"\n")
	file.close()

	file = open("bashrun.sh", "w")
	file.write("#!/bin/bash\n")
	if (seed == 1): file.write("rm -r "+folder+"\n")
	file.write("mkdir -p "+folder+ str(seed) +"\n")
	
	file.write("cp "+ compiled_fname+ " "+ folder + str(seed) + "\n")
	file.write("cp "+ job_script	+ " "+ folder + str(seed) + "\n")
	file.write("mv "+ input_i		+ " "+ folder + str(seed) +'/' + input_0 + "\n")
	
	for pfile in params:
		if os.path.exists(pfile): file.write("cp "+ pfile	+ " "+ folder + str(seed) + "\n")

	file.write("cd "+ folder + str(seed) + "\n")
	# file.write("python "+ compiled_fname + "\n")		
	file.write("qsub "+ job_script + "\n")
	file.write("cd .. \n")
	file.write("echo "+ "\"" + str(seed) + "\"" + " >> "+ dir_list +"\n") 
	file.write("cd .. \n")
	file.close()

	subprocess.call("(chmod +x bashrun.sh)", shell = True)
	subprocess.call("(./bashrun.sh)", shell = True)






