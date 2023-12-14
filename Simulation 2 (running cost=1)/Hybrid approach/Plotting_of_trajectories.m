%% Plotting of trajectories, hybrid method
% Juan Parras, GAPS-UPM, April 2016
clear all; clc; close all;
%% Load data
load('PAPER_error_oo_analytical_data_rc1_hybrid.mat');

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
tf=[];
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
                                param=struct('xe0',xe0(i1),'ye0',ye0(i2),'xp0',xp0(i3),'yp0',yp0(i4),'up0',up0(i5),'vp0',vp0(i6),'ue0',ue0(i7),'ve0',ve0(i8),'wp',wp,'Fp',Fp,'we',we,'Fe',Fe,'l',l,'ke',ke,'kp',kp,'D',D,'tf',tf);
                                s5_op=data_out(i1,i2,i3,i4,i5,i6,i7,i8,3);
                                s5_an=data_out(i1,i2,i3,i4,i5,i6,i7,i8,1);
                                tf_op=data_out(i1,i2,i3,i4,i5,i6,i7,i8,4);
                                tf_an=data_out(i1,i2,i3,i4,i5,i6,i7,i8,2);
                                param.tf=tf_op;
                                figure(1);
                                subplot(2,1,1);
                                if sum([isnan(s5_an),isnan(tf_an)])==0
                                    plot_traj_param_cap([s5_an tf_an],n_points);
                                else
                                    cla(subplot(2,1,1))
                                end
                                grid on;
                                subplot(2,1,2);
                                if sum([isnan(s5_op),isnan(tf_op)])==0
                                    plot_traj_param_cap([s5_op tf_op],n_points);
                                else
                                    cla(subplot(2,1,2))
                                end
                                grid on;
                                pause();
                            end
                        end
                    end
                end
            end
        end
    end
end
% %% Do the plotting
% 
% % Obtain grid points where solution was obtained
% points_correct_op=(out(:,:,:,:,:,:,:,:,4)>=-er_th);
% points_correct_op=points_correct_op(:); %Points well solved via hybrid method
% points_correct_an=(out(:,:,:,:,:,:,:,:,1)>0);
% points_correct_an=points_correct_an(:); %Points well solved via analytical method
% 
% s5_op=data_out(:,:,:,:,:,:,:,:,3);
% s5_an=data_out(:,:,:,:,:,:,:,:,1);
% tf_op=data_out(:,:,:,:,:,:,:,:,4);
% tf_an=data_out(:,:,:,:,:,:,:,:,2);
% s5_an=s5_an(:);
% s5_op=s5_op(:);
% tf_an=tf_an(:);
% tf_op=tf_op(:);
% 
% % Loop 
% 
% for i=1:length(parame)
%     %if xor(points_correct_an(i),points_correct_op(i)) % Enter loop if they differ
%         parame(i).tf=tf_op(i);
%         param=parame(i);
%         figure(1);
%         subplot(2,1,1);
%         plot_traj_param_cap([s5_an(i) tf_an(i)],n_points);
%         subplot(2,1,2);
%         plot_traj_param_cap([s5_op(i) tf_op(i)],n_points);
%         pause();
%     %end
% end
        