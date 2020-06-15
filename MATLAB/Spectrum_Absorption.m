%{
Created on Fri Oct 20 19:24:24 2017
@author: Lihua Xu  (libby)
%}

clc
close all %%Close all the plot, command "close all"----------------------%%
clear all
clc

%Scalling_factor=1000000;

%%load 9 files data to the newly made matrix-----------------------------%%
%%-----------------------------------------------------------------------%%
format long
load muxx.dat
length_step=length(muxx);
load muyx.dat
load muzx.dat
load muxy.dat
load muyy.dat
load muzy.dat
load muxz.dat
load muyz.dat
load muzz.dat
MuXx=muxx;
MuXy=muyx;
MuXz=muzx;
MuYx=muxy;
MuYy=muyy;
MuYz=muzy;
MuZx=muxz;
MuZy=muyz;
MuZz=muzz;

%%Shift the first value to 0 for each file (total 9)---------------------%%
%%-----------------------------------------------------------------------%%
Fisrt_value_Xx=MuXx(1,2);
Fisrt_value_Xy=MuXy(1,2);
Fisrt_value_Xz=MuXz(1,2);
Fisrt_value_Yx=MuYx(1,2);
Fisrt_value_Yy=MuYy(1,2);
Fisrt_value_Yz=MuYz(1,2);
Fisrt_value_Zx=MuZx(1,2);
Fisrt_value_Zy=MuZy(1,2);
Fisrt_value_Zz=MuZz(1,2);

MuXx(:,2)=MuXx(:,2)-Fisrt_value_Xx;
MuXy(:,2)=MuXy(:,2)-Fisrt_value_Xy;
MuXz(:,2)=MuXz(:,2)-Fisrt_value_Xz;
MuYx(:,2)=MuYx(:,2)-Fisrt_value_Yx;
MuYy(:,2)=MuYy(:,2)-Fisrt_value_Yy;
MuYz(:,2)=MuYz(:,2)-Fisrt_value_Yz;
MuZx(:,2)=MuZx(:,2)-Fisrt_value_Zx;
MuZy(:,2)=MuZy(:,2)-Fisrt_value_Zy;
MuZz(:,2)=MuZz(:,2)-Fisrt_value_Zz;

%%Average xx yy zz values from the files---------------------------------%%
%%-----------------------------------------------------------------------%%
Average_xxyyzz(:,1)=MuXx(:,1);
Average_xxyyzz(:,2)=(MuXx(:,2)+MuYy(:,2)+MuZz(:,2))./3;
figure(1)
set(gcf,'unit','normalized','position',[0.2,0.2,0.36,0.4])
plot(Average_xxyyzz(:,1),Average_xxyyzz(:,2),'-','Linewidth',1)
set(gca,'FontWeight','bold','Linewidth',1.3,'ticklength',2*get(gca,'ticklength'))
%axis auto
xlabel('Time (fs)','FontWeight','bold')
ylabel('Dipole moment','FontWeight','bold')
title('Before applying damping function','Fontsize',13)
%print('Before_applying_damping_function','-dpng','-r600');

%%Applying damping factor------------------------------------------------%%
%%-----------------------------------------------------------------------%%
dt=0.2; %atomic points
damp=zeros(length_step,1);
pdamp=0.008;
i=1:length_step;
j=i-1;
damp(i)=exp(-j.*dt*pdamp);
Average_xxyyzz(:,2)=Average_xxyyzz(:,2).*damp;
figure(2)
set(gcf,'unit','normalized','position',[0.2,0.2,0.36,0.4])
plot(Average_xxyyzz(:,1),Average_xxyyzz(:,2))
set(gca,'FontWeight','bold','Linewidth',1.3,'ticklength',2*get(gca,'ticklength'))
xlabel('Time (fs)','FontWeight','bold')
ylabel('Dipole moment','FontWeight','bold')
title('After applying damping function','Fontsize',13)
%print('After_applying_damping_function','-dpng','-r600');
%%You could add more zeros (8000 or more lines) at the end of the
%%Average_xxyyzz, This could help to increase the accuracy of the FTT
%%following.

%%Fourier transform------------------------------------------------------%%
%%-----------------------------------------------------------------------%%
%%Adding more data points into the data set (value equal to 0)-----------%%
temp=zeros(62013,2);
M=length(temp);
number_temp=length_step+1:length_step+M;
N=length_step+M;
temp(:,1)=(number_temp-1)*(MuXx(2,1)-MuXx(1,1));
Average_xxyyzz=[Average_xxyyzz;temp];
Average_xxyyzz_fft = fft(Average_xxyyzz(:,2),N);
imaginary_part=imag(Average_xxyyzz_fft);
imaginary_part=imaginary_part*(-1);
%length(imaginary_part)
%Finding what happened to the imaginary part
figure(3)
set(gcf,'unit','normalized','position',[0.2,0.2,0.36,0.4])
plot(1:length_step+M,imaginary_part)
%plot(1:length_step+M,imaginary_part.*Scalling_factor)
set(gca,'FontWeight','bold','Linewidth',1.3,'ticklength',2*get(gca,'ticklength'))
xlabel('Number','FontWeight','bold')
ylabel('Imaginary part','FontWeight','bold')
title('Temporary plot for imaginary part','Fontsize',13)
%print('Temporary_plot_for_imaginary_part','-dpng','-r600');
%%After FFT, only consider the first 5000 points, which is enough to get
%%the following results

%%Unit transfer----------------------------------------------------------%%
%%HZ---------------------------------------------------------------------%%
f=(0:(N-1))./(N*0.2); % Frequency atomic unit
Energy=f.*2*pi; % Energy atomic unit
Energy_ev=Energy.*27.2113961317875; % Energy eV unit
%length (Energy_ev)
figure(4)
set(gcf,'unit','normalized','position',[0.2,0.2,0.36,0.4])
plot(Energy_ev,imaginary_part,'Linewidth',2)
%plot(Energy_ev,imaginary_part.*Scalling_factor,'Linewidth',2)
set(gca,'FontWeight','bold','Linewidth',1.3,'ticklength',2*get(gca,'ticklength'))
xlabel('Energy (eV)','FontWeight','bold')
ylabel('Absorbance (arb.unit)','FontWeight','bold')
title('Absorbtance Spectrum','Fontsize',13)
axis([0 10 0 1.5])
%print('Absorbtance_Spectrum_Energy_Ev_Intensity','-dpng','-r600');

%%Unit transfer----------------------------------------------------------%%
%%Getting the Wavelength vs Energy------Experiment needed----------------%%
Wavelength=1239.8./Energy_ev;
figure(5)
set(gcf,'unit','normalized','position',[0.2,0.2,0.36,0.4])
plot(Wavelength,imaginary_part,'Linewidth',2)
%plot(Wavelength,imaginary_part.*Scalling_factor,'Linewidth',2)
set(gca,'FontWeight','bold','Linewidth',1.3,'ticklength',2*get(gca,'ticklength'))
xlabel('Wavelength (nm)','FontWeight','bold')
ylabel('Absorbance (arb.unit)','FontWeight','bold')
title('Absorbance Spectrum','Fontsize',13)
axis([100 1000 0 1.5])
%print('Absorbance_Spectrum_Wavelength_Intensity','-dpng','-r600');
%For Absorptance vs Energy(eV), We are interested in the energy range
%0-40eV.
%For Absorptance vs Wavelength (nm), We are interested in the wavelength
%range 100 to 1000 nm.

%%Applying the damping factor and FFT for all 9 cases--------------------%%
%%In order to find the polarizability vector-----------------------------%%
%%-----------------------------------------------------------------------%%
MuXx(:,2)=MuXx(:,2).*damp;
MuXy(:,2)=MuXy(:,2).*damp;
MuXz(:,2)=MuXz(:,2).*damp;
MuYx(:,2)=MuYx(:,2).*damp;
MuYy(:,2)=MuYy(:,2).*damp;
MuYz(:,2)=MuYz(:,2).*damp;
MuZx(:,2)=MuZx(:,2).*damp;
MuZy(:,2)=MuZy(:,2).*damp;
MuZz(:,2)=MuZz(:,2).*damp;

%%Applying the Fourier transform and FFT for all 9 cases-----------------%%
%%In order to find the polarizability vector-----------------------------%%
%%-----------------------------------------------------------------------%%
MuXx=[MuXx;temp];
MuXx_fft=fft(MuXx(:,2),N);
MuXx_imaginary_part=imag(MuXx_fft);
MuXx_imaginary_part=MuXx_imaginary_part*(-1);

MuXy=[MuXy;temp];
MuXy_fft=fft(MuXy(:,2),N);
MuXy_imaginary_part=imag(MuXy_fft);
MuXy_imaginary_part=MuXy_imaginary_part*(-1);

MuXz=[MuXz;temp];
MuXz_fft=fft(MuXz(:,2),N);
MuXz_imaginary_part=imag(MuXz_fft);
MuXz_imaginary_part=MuXz_imaginary_part*(-1);

MuYx=[MuYx;temp];
MuYx_fft=fft(MuYx(:,2),N);
MuYx_imaginary_part=imag(MuYx_fft);
MuYx_imaginary_part=MuYx_imaginary_part*(-1);

MuYy=[MuYy;temp];
MuYy_fft=fft(MuYy(:,2),N);
MuYy_imaginary_part=imag(MuYy_fft);
MuYy_imaginary_part=MuYy_imaginary_part*(-1);

MuYz=[MuYz;temp];
MuYz_fft=fft(MuYz(:,2),N);
MuYz_imaginary_part=imag(MuYz_fft);
MuYz_imaginary_part=MuYz_imaginary_part*(-1);

MuZx=[MuZx;temp];
MuZx_fft=fft(MuZx(:,2),N);
MuZx_imaginary_part=imag(MuZx_fft);
MuZx_imaginary_part=MuZx_imaginary_part*(-1);

MuZy=[MuZy;temp];
MuZy_fft=fft(MuZy(:,2),N);
MuZy_imaginary_part=imag(MuZy_fft);
MuZy_imaginary_part=MuZy_imaginary_part*(-1);

MuZz=[MuZz;temp];
MuZz_fft=fft(MuZz(:,2),N);
MuZz_imaginary_part=imag(MuZz_fft);
MuZz_imaginary_part=MuZz_imaginary_part*(-1);

%%Extract peak related to that in averaged value for these 9 cases-------%%
%%In order to find the polarizability vector-----------------------------%%
%%-----------------------------------------------------------------------%%
%%Find indexes for containning Energy_ev between 0 to 10-----------------%%
%%Dealing with the first 5000 points-------------------------------------%%
Energy_ev_temp=Energy_ev(1:5000);
indexes=find(Energy_ev_temp>=0&Energy_ev_temp<=10);
%%Energy we are interested are always within this range------------------%%
small_value=min(indexes);
big_value=max(indexes);
%%Find indexes for the peak for averged data among xx yy zz--------------%%
%max(imaginary_part(small_value:big_value))
max_peak_index_Ave=find(imaginary_part==max(imaginary_part(small_value:big_value)));
%%Disp the Energy_ev for the max_peak_index_Ave
disp('The Energy (eV) for the max peak during this period is')
Energy_ev(max_peak_index_Ave)
%%The range to find the largest peak (intensity)


%%This is using the lower triangle as the base
matrix_nine=zeros(3,3);
matrix_nine(2,1)=MuYx_imaginary_part(max_peak_index_Ave);
matrix_nine(3,1)=MuZx_imaginary_part(max_peak_index_Ave);
matrix_nine(3,2)=MuZy_imaginary_part(max_peak_index_Ave);
matrix_nine(1,1)=MuXx_imaginary_part(max_peak_index_Ave);
matrix_nine(2,2)=MuYy_imaginary_part(max_peak_index_Ave);
matrix_nine(3,3)=MuZz_imaginary_part(max_peak_index_Ave);
matrix_nine(1,2)=matrix_nine(2,1);
matrix_nine(2,3)=matrix_nine(3,2);
matrix_nine(1,3)=matrix_nine(3,1);
%matrix_nine

%{ 
%%This is using the upper triangle as the base
matrix_nine(1,1)=MuXx_imaginary_part(max_peak_index_Ave);
matrix_nine(1,2)=MuXy_imaginary_part(max_peak_index_Ave);
matrix_nine(1,3)=MuXz_imaginary_part(max_peak_index_Ave);
matrix_nine(2,1)=matrix_nine(1,2);
%matrix_nine(2,1)=MuYx_imaginary_part(max_peak_index_Ave);
matrix_nine(2,2)=MuYy_imaginary_part(max_peak_index_Ave);
matrix_nine(2,3)=MuYz_imaginary_part(max_peak_index_Ave);
matrix_nine(3,1)=matrix_nine(1,3);
%matrix_nine(3,1)=MuZx_imaginary_part(max_peak_index_Ave);
matrix_nine(3,2)=matrix_nine(2,3);
%matrix_nine(3,2)=MuZy_imaginary_part(max_peak_index_Ave);
matrix_nine(3,3)=MuZz_imaginary_part(max_peak_index_Ave);
%}

%%Getting the eigenvalues and eigenvectors for the above matrix----------%%
%%In order to find the polarizability vector-----------------------------%%
%%-----------------------------------------------------------------------%%
[Eigenvectors,Eigenvalues] = eig(matrix_nine);
%%Finding the maxmum eigenvalue
value_in_matrix=max(max(Eigenvalues));
[row col]=find(value_in_matrix==Eigenvalues);
disp('1-x, 2-y,3-z,so the largest number for the eigenvalues is')
row
disp('The conrresponding eigenvectors is')
Eigenvectors(:,row)

%close all 
%%clear all
%%clc
