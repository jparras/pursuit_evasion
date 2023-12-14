%% Plotting of trajectories, hybrid method
% Juan Parras, GAPS-UPM, April 2016
clear all; clc; close all;
%% Load data
load('PAPER_error_oo_data_rcr.mat');

%% Parameters initialization
xe0=1:5:11;
ye0=1:5:11;
xp0=-10:5:0;
yp0=-10:5:0;
up0=-1;
vp0=-1;
ue0=1;
ve0=1;
wp=2;
Fp=1;
we=1;
Fe=1;
l=1;
ke=Fe/we;
kp=Fp/wp;
D=100; %Maximum abs value of the coordinates
tf_max=4*D/(wp-we); %Sufficiently high
Pj=1.11;
Pc=1;
N=100;
A=N*(log2(1+Pc/Pj)+(Pj*(0.0069*D*(1+Pj/Pc)/sqrt(Pj/Pc)+14.4070)/Pc)/(2*D^2*(1+Pj/Pc)^2*log(2)));
B=N*(Pj*log(0.1824*D*(1+Pj/Pc)/sqrt(Pj/Pc)+0.4823)/Pc)/(2*D^2*(1+Pj/Pc)^2*log(2));
dv=2*max(wp,we);
denorm_vector=[2*pi tf_max dv];
er_th=0.9;
global param; % Declared as global, to be used by my_cost_fun
n_points=1e3;

for i1=1:length(xe0)
    for i2=1:length(ye0)
        for i3=1:length(xp0)
            for i4=1:length(yp0)
                for i5=1:length(up0)
                    for i6=1:length(vp0)
                        for i7=1:length(ue0)
                            for i8=1:length(ve0)
                                param=struct('xe0',xe0(i1),'ye0',ye0(i2),'xp0',xp0(i3),'yp0',yp0(i4),'up0',up0(i5),'vp0',vp0(i6),'ue0',ue0(i7),'ve0',ve0(i8),'wp',wp,'Fp',Fp,'we',we,'Fe',Fe,'l',l,'ke',ke,'kp',kp,'D',D,'tf_max',tf_max,'A',A,'B',B,'dv',dv);
                                s5=data_out(i1,i2,i3,i4,i5,i6,i7,i8,1);
                                tf=data_out(i1,i2,i3,i4,i5,i6,i7,i8,2);
                                dv=data_out(i1,i2,i3,i4,i5,i6,i7,i8,3);
                                figure(1);
                                if out(i1,i2,i3,i4,i5,i6,i7,i8,1)>=-er_th
                                    plot_traj_param_rcr([s5 tf dv],n_points);
                                else
                                    cla(figure(1))
                                end
                                % To see all trajectories, get the pause
                                % out of the if
                                if param.xe0==6&&param.ye0==1&&param.xp0==-10&&param.yp0==0
                                    pause();
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end