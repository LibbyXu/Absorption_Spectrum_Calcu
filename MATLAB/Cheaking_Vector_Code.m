%{
Created on Fri Oct 20 19:24:24 2017
@author: Lihua Xu  (libby)
%}

clc
clear all
close all
clc

%%Recheck the peak information produced by the former Matlab Codes-------%%
%%-----------------------------------------------------------------------%%
format long
load mux.dat
length_step=length(mux);
load muy.dat
load muz.dat
MuX=mux;
MuY=muy;
MuZ=muz;

%%Shift the first value to 0 for each file (total 3)---------------------%%
%%-----------------------------------------------------------------------%%
Fisrt_value_X=MuX(1,2);
Fisrt_value_Y=MuY(1,2);
Fisrt_value_Z=MuY(1,2);

%%Norming the dipole moment and then get the plot------------------------%%
%%-----------------------------------------------------------------------%%
MuX(:,2)=MuX(:,2)-Fisrt_value_X;
MuY(:,2)=MuY(:,2)-Fisrt_value_Y;
MuZ(:,2)=MuZ(:,2)-Fisrt_value_Z;

Mu(:,2)=sqrt(MuX(:,2).^2+MuY(:,2).^2+MuZ(:,2).^2);
Mu(:,1)=MuX(:,1);

figure(6)
set(gcf,'unit','normalized','position',[0.2,0.2,0.36,0.4])
plot(MuX(:,1),MuX(:,2),'-b','Linewidth',1.2)
set(gca,'FontWeight','bold','Linewidth',1.3,'ticklength',2*get(gca,'ticklength'))
hold on
plot(MuY(:,1),MuY(:,2),'-r','Linewidth',1.2)
plot(MuZ(:,1),MuZ(:,2),'-g','Linewidth',1.2)
xlabel('Time (fs)','FontWeight','bold')
ylabel('Dipole moment','FontWeight','bold')
%title('','Fontsize',13)
legend({'Blue-MuX','Red-MuY','Green-MuZ'},'Fontsize',10,'Location','northwest')
hold off

figure(10)
set(gcf,'unit','normalized','position',[0.2,0.2,0.36,0.4])
plot(Mu(1:10000,1),Mu(1:10000,2),'-','Linewidth',1)
%%The plot should have a linearity increasing tendency-------------------%%
set(gca,'FontWeight','bold','Linewidth',1.3,'ticklength',2*get(gca,'ticklength'))
%axis auto
xlabel('Time (fs)','FontWeight','bold')
ylabel('Norm Dipole moment','FontWeight','bold')
title('Cheak the peak produced by the former Matlab code','Fontsize',13)
%print('Cheak_the_peak_produced_by_the_former_Matlab_code','-dpng','-r600');

%{
figure(7)
set(gcf,'unit','normalized','position',[0.2,0.2,0.36,0.4])
plot(MuX(1:10000,1),MuX(1:10000,2),'-','Linewidth',1)
%%The plot should have a linearity increasing tendency-------------------%%
set(gca,'FontWeight','bold','Linewidth',1.3,'ticklength',2*get(gca,'ticklength'))
%axis auto
xlabel('Time (fs)','FontWeight','bold')
ylabel('X Direction Dipole moment','FontWeight','bold')
title('X','Fontsize',13)
%print('Cheak_the_peak_produced_by_the_former_Matlab_code','-dpng','-r600');

figure(8)
set(gcf,'unit','normalized','position',[0.2,0.2,0.36,0.4])
plot(MuY(1:10000,1),MuY(1:10000,2),'-','Linewidth',1)
%%The plot should have a linearity increasing tendency-------------------%%
set(gca,'FontWeight','bold','Linewidth',1.3,'ticklength',2*get(gca,'ticklength'))
%axis auto
xlabel('Time (fs)','FontWeight','bold')
ylabel('Y Direction Dipole moment','FontWeight','bold')
title('Y','Fontsize',13)
%print('Cheak_the_peak_produced_by_the_former_Matlab_code','-dpng','-r600');

figure(9)
set(gcf,'unit','normalized','position',[0.2,0.2,0.36,0.4])
plot(MuZ(1:10000,1),MuZ(1:10000,2),'-','Linewidth',1)
%%The plot should have a linearity increasing tendency-------------------%%
set(gca,'FontWeight','bold','Linewidth',1.3,'ticklength',2*get(gca,'ticklength'))
%axis auto
xlabel('Time (fs)','FontWeight','bold')
ylabel('Z Direction Dipole moment','FontWeight','bold')
title('Z','Fontsize',13)
%print('Cheak_the_peak_produced_by_the_former_Matlab_code','-dpng','-r600');
%}