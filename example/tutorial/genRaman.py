import random
import numpy as np

length = 100000
oscillators = 2
xx = 1 #20.0
xy = 0.0
xz = 0.0
yy = 1 #4.0
yz = 0.0
zz = 1.0

#raman=open("RamanTensor.txt","w")
file_raman=open("RamanTensor.bin","wb")
for i in range(int(length)):
   step=np.array([0],'float32')
   step.tofile(file_raman)
   R4bin=np.array([xx,xx,xy,xy,xz,xz,yy,yy,yz,yz,zz,zz])
   Rf=np.array(R4bin,'float32')
   Rf.tofile(file_raman)
  # rewrite that is prints for singles
#  raman.write(str(i) + " "+ str(xx) + " " + str(xx) + " " + str(xy) + " "+ str(xy) + " "  + str(xz) + " " + str(xz) + " " + str(yy) + " " + str(yy) + " " + str(yz) + " " + str(yz) + " "+ str(zz) + " "  + str(zz) + "\n")

file_raman.close()
