import os
import numpy as np

folder = "./run1/"
dir_list_file = "./run1/pfolder.txt"
file = "/trajStates.dat"

momentums = np.loadtxt(dir_list_file,dtype=str) #os.listdir(folder)
for pval in momentums:
  p, s, t = np.loadtxt(folder+pval+file, unpack=True)
  print(t)
  ttotal = t[0]+t[1]+t[3]+t[4]
  with open(folder+"lu.dat", 'a') as f0 : print(p[0],s[0],t[0]/ttotal, file=f0)
  with open(folder+"ld.dat", 'a') as f1 : print(p[1],s[1],t[1]/ttotal, file=f1)
  with open(folder+"zero.dat",'a') as f2: print(p[2],s[2],t[2]/ttotal, file=f2)
  with open(folder+"ru.dat", 'a') as f3 : print(p[3],s[3],t[3]/ttotal, file=f3)
  with open(folder+"rd.dat", 'a') as f4 : print(p[4],s[4],t[4]/ttotal, file=f4)

#   np.savetxt('lu.dat', np.c_[ p[0],s[0],t[0] ], a)
#   np.savetxt('ld.dat', np.c_[ p[1],s[1],t[1] ], a)
#   np.savetxt('rd.dat', np.c_[ p[3],s[3],t[3] ], a)
#   np.savetxt('ru.dat', np.c_[ p[4],s[4],t[4] ], a)

