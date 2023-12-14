%% Comparison of controls of capacity and pursuit-evasion game
% Juan Parras, GAPS-UPM, April 2016
clear all; clc; close all;
%% Load data from simulation 2
load('PAPER_error_oo_analytical_data_rc1_hybrid.mat');
out_rc1=out;
data_out_rc1=data_out;
%% Load data from simulation 3
load('PAPER_error_oo_data_rcr_dv_appr.mat');
out_rcr=out;
data_out_rcr=data_out;
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
Pj=1.11;
Pc=1;
N=100;
A=N*(log2(1+Pc/Pj)+(Pj*(0.0069*D*(1+Pj/Pc)/sqrt(Pj/Pc)+14.4070)/Pc)/(2*D^2*(1+Pj/Pc)^2*log(2)));
B=N*(Pj*log(0.1824*D*(1+Pj/Pc)/sqrt(Pj/Pc)+0.4823)/Pc)/(2*D^2*(1+Pj/Pc)^2*log(2));
dv=2*max(wp,we);
global param; % Declared as global, to be used by my_cost_fun

n_points=1e3;

er_th=0.9;

pl=0; %To plot in each iteration the trajectories obtained
n_p_t=1e2; %Points in plot trajectories

%% Loop
i=0;
if pl
    figure(1);
end
error_p=[];
error_e=[];
for i1=1:length(xe0)
    for i2=1:length(ye0)
        for i3=1:length(xp0)
            for i4=1:length(yp0)
                for i5=1:length(up0)
                    for i6=1:length(vp0)
                        for i7=1:length(ue0)
                            for i8=1:length(ve0)
                                % Check that RCR finds solution
                                if out_rcr(i1,i2,i3,i4,i5,i6,i7,i8,1)>=-er_th
                                    i=i+1;
                                    param=struct('xe0',xe0(i1),'ye0',ye0(i2),'xp0',xp0(i3),'yp0',yp0(i4),'up0',up0(i5),'vp0',vp0(i6),'ue0',ue0(i7),'ve0',ve0(i8),'wp',wp,'Fp',Fp,'we',we,'Fe',Fe,'l',l,'ke',ke,'kp',kp,'D',D,'tf',tf,'tf_max',tf_max,'A',A,'B',B,'dv',dv);
                                    % Control with running cost=1 (hybrid approach)
                                    control_rc1(i)=data_out_rc1(i1,i2,i3,i4,i5,i6,i7,i8,3); %control is constant!!
                                    tf_rc1(i)=data_out_rc1(i1,i2,i3,i4,i5,i6,i7,i8,4);
                                    
                                    if pl
                                        plot_traj_param_cap([control_rc1(i) tf_rc1(i)],n_p_t);
                                    end
                                    % Controls with running cost=r 
                                    s5=data_out_rcr(i1,i2,i3,i4,i5,i6,i7,i8,1);
                                    tf=data_out_rcr(i1,i2,i3,i4,i5,i6,i7,i8,2);
                                    delta_v=wp-we;
                                    
                                    [xs,t]=obtain_xs(s5,tf,delta_v,n_points);

                                    xp=xs(:,1);
                                    up=xs(:,2);
                                    yp=xs(:,3);
                                    vp=xs(:,4);
                                    xe=xs(:,5);
                                    ue=xs(:,6);
                                    ye=xs(:,7);
                                    ve=xs(:,8);
                                    
                                    if pl
                                        plot_traj_param_rcr([s5 tf delta_v],n_p_t);
                                        pause();
                                    end

                                    Bphi=eps+2*B*(exp(kp*(tf-t)).*(kp*(tf-t)-1)+1).*(xp-xe)+kp*((A+B*((yp-ye).^2+(xp-xe).^2))./delta_v).*(exp(kp*(tf-t))-1).*sin(s5);
                                    Aphi=eps+2*B*(exp(kp*(tf-t)).*(kp*(tf-t)-1)+1).*(yp-ye)+kp*((A+B*((yp-ye).^2+(xp-xe).^2))./delta_v).*(exp(kp*(tf-t))-1).*cos(s5);
                                    Bpsi=eps+2*B*(exp(ke*(tf-t)).*(ke*(tf-t)-1)+1).*(xp-xe)+ke*((A+B*((yp-ye).^2+(xp-xe).^2))./delta_v).*(exp(ke*(tf-t))-1).*sin(s5);
                                    Apsi=eps+2*B*(exp(ke*(tf-t)).*(ke*(tf-t)-1)+1).*(yp-ye)+ke*((A+B*((yp-ye).^2+(xp-xe).^2))./delta_v).*(exp(ke*(tf-t))-1).*cos(s5);

                                    cos_phi=Aphi./sqrt(eps+Aphi.^2+Bphi.^2);
                                    sin_phi=Bphi./sqrt(eps+Aphi.^2+Bphi.^2);
                                    cos_psi=Apsi./sqrt(eps+Apsi.^2+Bpsi.^2);
                                    sin_psi=Bpsi./sqrt(eps+Apsi.^2+Bpsi.^2);
                                    
                                    ang_phi=atan2(sin_phi,cos_phi);
                                    ang_psi=atan2(sin_psi,cos_psi);
                                    
                                    err_abs_pursuer=abs(mod(ang_phi, 2*pi)-mod(control_rc1(i), 2*pi));
                                    err_abs_evader=abs(mod(ang_psi, 2*pi)-mod(control_rc1(i), 2*pi));
                                    
                                    err_rel_pursuer=err_abs_pursuer./(eps+abs(mod(ang_phi, 2*pi)));
                                    err_rel_evader=err_abs_evader./(eps+abs(mod(ang_psi, 2*pi)));
                                    
                                    error_p=[error_p; err_rel_pursuer];
                                    error_e=[error_e; err_rel_evader];
                                    
                                    e_rel_p_mean(i)=mean(err_rel_pursuer);
                                    e_rel_e_mean(i)=mean(err_rel_evader);
                                                                        
                                    e_rel_p_median(i)=median(err_rel_pursuer);
                                    e_rel_e_median(i)=median(err_rel_evader);
                                    
                                    e_rel_p_var(i)=var(err_rel_pursuer);
                                    e_rel_e_var(i)=var(err_rel_evader);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Plot

figure(2);
subplot(3,1,1);
plot(1:i,e_rel_p_mean*100,'r',1:i,e_rel_e_mean*100,'b');
xlabel('Grid point');
ylabel('Mean relative error (%)');
grid on;
subplot(3,1,2);
plot(1:i,e_rel_p_median*100,'r',1:i,e_rel_e_median*100,'b');
xlabel('Grid point');
ylabel('Median relative error (%)');
grid on;
subplot(3,1,3);
plot(1:i,100*sqrt(e_rel_p_var),'r',1:i,100*sqrt(e_rel_e_var),'b');
xlabel('Grid point');
ylabel('Relative error standard deviation (%)');
grid on;

display(['Pursuer total error metrics: mean(%) = ' num2str(100*mean(error_p)) ' ; median(%) = ' num2str(100*median(error_p)) ' ; standard deviation(%) = ' num2str(100*sqrt(var(error_p)))]);
display(['Evader total error metrics: mean(%) = ' num2str(100*mean(error_e)) ' ; median(%) = ' num2str(100*median(error_e)) ' ; standard deviation(%) = ' num2str(100*sqrt(var(error_e)))]);
figure(3);
subplot(1,2,1);
histogram(100*error_e,100,'Normalization','Probability');
title('Evader error hist (%)');
subplot(1,2,2);
histogram(100*error_p,100,'Normalization','Probability');
title('Pursuer error hist (%)');
figure(4);
subplot(1,2,1);
cdfplot(100*error_e);
subplot(1,2,2);
cdfplot(100*error_p);
%% Save CDF values
[y_e_delta, x_e_delta]=ecdf(100*error_e);
[y_p_delta, x_p_delta]=ecdf(100*error_p);
save('CDF_delta','y_e_delta','x_e_delta','y_p_delta','x_p_delta');
