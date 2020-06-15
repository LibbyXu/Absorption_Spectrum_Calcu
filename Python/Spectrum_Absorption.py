"""
Created on Thu Oct 19 16:21:33 2017
@author: Lihua Xu (libby)
"""

#Scalling_factor=1000000;

##load 9 files data to the newly made matrix
import numpy as np
MuXx = np.genfromtxt('muxx.dat',delimiter='')
MuXy = np.genfromtxt('muyx.dat',delimiter='')
MuXz = np.genfromtxt('muzx.dat',delimiter='')
MuYx = np.genfromtxt('muxy.dat',delimiter='')
MuYy = np.genfromtxt('muyy.dat',delimiter='')
MuYz = np.genfromtxt('muzy.dat',delimiter='')
MuZx = np.genfromtxt('muxz.dat',delimiter='')
MuZy = np.genfromtxt('muyz.dat',delimiter='')
MuZz = np.genfromtxt('muzz.dat',delimiter='')

#print(MuXx)
#print(MuXy)
#print(MuXz)
#print(MuYx)
#print(MuYy)
#print(MuYz)
#print(MuZx)
#print(MuZy)
#print(MuZz)

Fisrt_value_Xx=MuXx[0][1]
Fisrt_value_Xy=MuXy[0][1]
Fisrt_value_Xz=MuXz[0][1]
Fisrt_value_Yx=MuYx[0][1]
Fisrt_value_Yy=MuYy[0][1]
Fisrt_value_Yz=MuYz[0][1]
Fisrt_value_Zx=MuZx[0][1]
Fisrt_value_Zy=MuZy[0][1]
Fisrt_value_Zz=MuZz[0][1]

##Shift the second column in order to make the first value be 0

MuXx_mu=MuXx[:,1]-Fisrt_value_Xx
MuXy_mu=MuXy[:,1]-Fisrt_value_Xy
MuXz_mu=MuXz[:,1]-Fisrt_value_Xz
MuYx_mu=MuYx[:,1]-Fisrt_value_Yx
MuYy_mu=MuYy[:,1]-Fisrt_value_Yy
MuYz_mu=MuYz[:,1]-Fisrt_value_Yz
MuZx_mu=MuZx[:,1]-Fisrt_value_Zx
MuZy_mu=MuZy[:,1]-Fisrt_value_Zy
MuZz_mu=MuZz[:,1]-Fisrt_value_Zz
length_step=len(MuZz_mu)
Mu_num=np.zeros((length_step,1))
for t in range(0,length_step):
    Mu_num[t]=MuZz[t][0]

#Average xx yy zz values from the files
Average_xxyyzz=(MuXx_mu+MuYy_mu+MuZz_mu)/3
import matplotlib.pylab as pl
pl.figure(1)
#pl.subplot(311)
pl.plot(Mu_num, Average_xxyyzz)
pl.title('Before applying damping function')# give plot a title
pl.xlabel('Time (fs)')# make axis labels
pl.ylabel('Dipole moment')
pl.xlim(0, 100)# set axis limits
#pl.ylim(0, 1)
pl.show()

#Applying damping factor
dt=0.2
pdamp=0.008
damp=np.zeros((length_step,1))
Average_xxyyzz_new=np.zeros((length_step,1))
for i in range(0,length_step):
    damp[i]=np.exp(-i*dt*pdamp)
for j in range(0,length_step):
    Average_xxyyzz_new[j]=Average_xxyyzz[j]*damp[j]
pl.figure(2)
#pl.subplot(312)
pl.plot(Mu_num, Average_xxyyzz_new)
pl.title('After applying damping function')# give plot a title
pl.xlabel('Time (fs)')# make axis labels
pl.ylabel('Dipole moment')
pl.xlim(0, 100)# set axis limits    
#pl.ylim(0, 1)
pl.show()

#Fourier transform
temp_mu=np.zeros(((length_step-1)*3,1))
temp_N=len(temp_mu)
Average_extended=np.vstack((Average_xxyyzz_new,temp_mu))
N=temp_N+length_step
temp_num=np.zeros(((length_step-1)*3,1))
for mm in range(length_step,N):
    temp_num[mm-length_step]=mm*(Mu_num[1]-Mu_num[0])
Mu_extended=np.vstack((Mu_num,temp_num))
##Do Fourier transform
from scipy.fftpack import fft
Average_extended_new=Average_extended[:,0]
Average_extended_fft=fft(Average_extended_new)
Average_imaginary = Average_extended_fft.imag
Average_imaginary = Average_imaginary*(-1)
order=np.zeros((N,1))
for ni in range(0,N):
    order[ni]=ni  
pl.figure(3)
#pl.subplot(313)
pl.plot(order, Average_imaginary)
pl.title('Temporary plot for imaginary part')# give plot a title
pl.xlabel('Number')# make axis labels
pl.ylabel('Imaginary part')
pl.xlim(0, N)# set axis limits 
pl.show()

#Unit transfer 
##HZ
import math
f=np.zeros((N,1))
Energy=np.zeros((N,1))
for Ut in range(0,N):
    f[Ut]=Ut/(dt*N)
    Energy[Ut]=f[Ut]*2*math.pi
Hartree_ev=27.2113961317875
Energy_ev=Energy*Hartree_ev
pl.figure(4)
pl.plot(Energy_ev, Average_imaginary)
pl.title('Absorbtance Spectrum')# give plot a title
pl.xlabel('Energy (eV)')# make axis labels
pl.ylabel('Absorbance (arb.unit)')
pl.xlim(0, 10)# set axis limits    
#pl.ylim(0, 1)
pl.show()

##Getting the Wavelength vs Energy  (Experiment needed)
Wavelength=1239.8/Energy_ev
pl.figure(5)
pl.plot(Wavelength, Average_imaginary)
pl.title('Absorbtance Spectrum')# give plot a title
pl.xlabel('Wavelength (nm)')# make axis labels
pl.ylabel('Absorbance (arb.unit)')
pl.xlim(100, 1000)# set axis limits    
#pl.ylim(0, 1)
pl.show()

#Applying the damping factor and FFT for all 9 cases
##In order to find the polarizability vector
MuXx_mu_damp=np.zeros((length_step,1))
MuXy_mu_damp=np.zeros((length_step,1))
MuXz_mu_damp=np.zeros((length_step,1))
MuYx_mu_damp=np.zeros((length_step,1))
MuYy_mu_damp=np.zeros((length_step,1))
MuYz_mu_damp=np.zeros((length_step,1))
MuZx_mu_damp=np.zeros((length_step,1))
MuZy_mu_damp=np.zeros((length_step,1))
MuZz_mu_damp=np.zeros((length_step,1))
for jj in range(0,length_step):
    MuXx_mu_damp[jj]=MuXx_mu[jj]*damp[jj]
    MuXy_mu_damp[jj]=MuXy_mu[jj]*damp[jj]
    MuXz_mu_damp[jj]=MuXz_mu[jj]*damp[jj]
    MuYx_mu_damp[jj]=MuYx_mu[jj]*damp[jj]
    MuYy_mu_damp[jj]=MuYy_mu[jj]*damp[jj]
    MuYz_mu_damp[jj]=MuYz_mu[jj]*damp[jj]
    MuZx_mu_damp[jj]=MuZx_mu[jj]*damp[jj]
    MuZy_mu_damp[jj]=MuZy_mu[jj]*damp[jj]
    MuZz_mu_damp[jj]=MuZz_mu[jj]*damp[jj]
    
##In order to find the polarizability vector
MuXx_extended=np.vstack((MuXx_mu_damp,temp_mu))
MuXy_extended=np.vstack((MuXy_mu_damp,temp_mu))
MuXz_extended=np.vstack((MuXz_mu_damp,temp_mu))
MuYx_extended=np.vstack((MuYx_mu_damp,temp_mu))
MuYy_extended=np.vstack((MuYy_mu_damp,temp_mu))
MuYz_extended=np.vstack((MuYz_mu_damp,temp_mu))
MuZx_extended=np.vstack((MuZx_mu_damp,temp_mu))
MuZy_extended=np.vstack((MuZy_mu_damp,temp_mu))
MuZz_extended=np.vstack((MuZz_mu_damp,temp_mu))

MuXx_extended_new=MuXx_extended[:,0]
MuXx_extended_fft=fft(MuXx_extended_new)
MuXx_imaginary = MuXx_extended_fft.imag
MuXx_imaginary = MuXx_imaginary*(-1)

MuXy_extended_new=MuXy_extended[:,0]
MuXy_extended_fft=fft(MuXy_extended_new)
MuXy_imaginary = MuXy_extended_fft.imag
MuXy_imaginary = MuXy_imaginary*(-1)

MuXz_extended_new=MuXz_extended[:,0]
MuXz_extended_fft=fft(MuXz_extended_new)
MuXz_imaginary = MuXz_extended_fft.imag
MuXz_imaginary = MuXz_imaginary*(-1)

MuYx_extended_new=MuYx_extended[:,0]
MuYx_extended_fft=fft(MuYx_extended_new)
MuYx_imaginary = MuYx_extended_fft.imag
MuYx_imaginary = MuYx_imaginary*(-1)

MuYy_extended_new=MuYy_extended[:,0]
MuYy_extended_fft=fft(MuYy_extended_new)
MuYy_imaginary = MuYy_extended_fft.imag
MuYy_imaginary = MuYy_imaginary*(-1)

MuYz_extended_new=MuYz_extended[:,0]
MuYz_extended_fft=fft(MuYz_extended_new)
MuYz_imaginary = MuYz_extended_fft.imag
MuYz_imaginary = MuYz_imaginary*(-1)

MuZx_extended_new=MuZx_extended[:,0]
MuZx_extended_fft=fft(MuZx_extended_new)
MuZx_imaginary = MuZx_extended_fft.imag
MuZx_imaginary = MuZx_imaginary*(-1)

MuZy_extended_new=MuZy_extended[:,0]
MuZy_extended_fft=fft(MuZy_extended_new)
MuZy_imaginary = MuZy_extended_fft.imag
MuZy_imaginary = MuZy_imaginary*(-1)

MuZz_extended_new=MuZz_extended[:,0]
MuZz_extended_fft=fft(MuZz_extended_new)
MuZz_imaginary = MuZz_extended_fft.imag
MuZz_imaginary = MuZz_imaginary*(-1)

#Extract peak related to that in averaged value for these 9 cases-
index_vector=np.zeros((5000,1))
iv=len(index_vector)
pp=0
for E_num in range(0,5000):
    if Energy_ev[E_num]>=0 and Energy_ev[E_num]<=10:
        index_vector[pp]=E_num
        pp=pp+1
import scipy
index_vector=scipy.delete(index_vector,range(pp,iv),0)
small_value=min(index_vector)
small_value=int(small_value[0])
big_value=max(index_vector)
big_value=int(big_value[0])
inteval=big_value-small_value+1
Average_imaginary_new=np.zeros((N,1))
for cc in range(0,N):
    Average_imaginary_new[cc]=Average_imaginary[cc]
existing_matrix_index=np.zeros((inteval,1))
existing_matrix_value=np.zeros((inteval,1))
for ff in range(0,inteval):
    existing_matrix_value[ff,0]=Average_imaginary_new[small_value+ff]
    existing_matrix_index[ff,0]=small_value+ff
peak=max(existing_matrix_value)
#peak=peak[0]
existing_matrix_value_list = existing_matrix_value.tolist()
peak_index=existing_matrix_index[existing_matrix_value_list.index(peak),0]
peak_index=int(peak_index)
print "The Energy (eV) for the max peak during this period is "
Peak_Energy=Energy_ev[peak_index][0]
print Peak_Energy

##This is using the lower triangle as the base
matrix_nine=np.zeros((3,3))
matrix_nine[0][0]=MuXx_imaginary[peak_index]
matrix_nine[1][0]=MuYx_imaginary[peak_index]
matrix_nine[2][0]=MuZx_imaginary[peak_index]
matrix_nine[2][1]=MuZy_imaginary[peak_index]
matrix_nine[1][1]=MuYy_imaginary[peak_index]
matrix_nine[2][2]=MuZz_imaginary[peak_index]
matrix_nine[0][1]=matrix_nine[1][0]
matrix_nine[1][2]=matrix_nine[2][1]
matrix_nine[0][2]=matrix_nine[2][0]
#print matrix_nine

#Getting the eigenvalues and eigenvectors for the above matrix
#In order to find the polarizability vector
from numpy import linalg as LA
w, v = LA.eig(np.array(matrix_nine))
Max_eigenvalue=max(w)
Eigenvalues = w.tolist()
Max_Eigvalue_index=Eigenvalues.index(Max_eigenvalue)
print "The maximum eigenvalue is "
print Eigenvalues[Max_Eigvalue_index]
print "The conrresponding eigenvectors is for x y and z direction are "
print v[:,Max_Eigvalue_index]
