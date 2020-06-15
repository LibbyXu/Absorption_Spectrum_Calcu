"""
Created on Fri Oct 20 19:24:24 2017
@author: Lihua Xu  (libby)
"""

##load 3 files data to the newly made matrix
##Recheck the peak information produced by the former Matlab Codes
import numpy as np
MuX = np.genfromtxt('mux.dat',delimiter='')
MuY = np.genfromtxt('muy.dat',delimiter='')
MuZ = np.genfromtxt('muz.dat',delimiter='')

##Shift the first value to 0 for each file (total 3)
Fisrt_value_X=MuX[0][1]
Fisrt_value_Y=MuY[0][1]
Fisrt_value_Z=MuZ[0][1]

##Norming the dipole moment and then get the plot
MuX_mu=MuX[:,1]-Fisrt_value_X
MuY_mu=MuY[:,1]-Fisrt_value_Y
MuZ_mu=MuZ[:,1]-Fisrt_value_Z
length_step=len(MuZ_mu)
Mu_num=np.zeros((length_step,1))
for t in range(0,length_step):
    Mu_num[t]=MuZ[t][0]

#MuX_mu_list = MuX_mu.tolist()
#MuY_mu_list = MuY_mu.tolist()
#MuZ_mu_list = MuZ_mu.tolist() 
XX=np.square(MuX_mu)
YY=np.square(MuY_mu)
ZZ=np.square(MuZ_mu)
norm_xyz=np.sqrt(XX+YY+ZZ)
#print XX
#print YY
#print ZZ

import matplotlib.pylab as pl
pl.figure(1)
pl.plot(Mu_num, norm_xyz)
pl.title('Cheak the peak produced by the former Matlab code')# give plot a title
pl.xlabel('Time (fs)')# make axis labels
pl.ylabel('Norm Dipole moment')
#pl.xlim(0, 100)# set axis limits
#pl.ylim(0, 1)
pl.show()

pl.figure(2)
label = ["MuX","MuY","MuZ"] 
pl.plot(Mu_num, MuX_mu,'r')
pl.plot(Mu_num, MuY_mu,'b')
pl.plot(Mu_num, MuZ_mu,'g')
pl.legend(label,loc='best')
pl.xlabel('Time (fs)')# make axis labels
pl.ylabel('Dipole moment')
#pl.xlim(0, 100)# set axis limits
#pl.ylim(0, 1)
pl.show()
