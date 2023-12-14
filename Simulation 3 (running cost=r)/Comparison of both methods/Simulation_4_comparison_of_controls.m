%% Comparison of controls of capacity game optimization and delta_v approaches
% Juan Parras, GAPS-UPM, April 2016
clear all; clc; close all;
%% Load data from simulation 3, optimization approach
load('PAPER_error_oo_data_rcr.mat');
out_rcr=out;
data_out_rcr=data_out;
%% Load data from simulation 3, delta_v approach
load('PAPER_error_oo_data_rcr_dv_appr.mat');
out_rcr2=out;
data_out_rcr2=data_out;
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

er_th=0.9;


%% Loop
i=0;
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
                                % Check that there are common solutions
                                if out_rcr(i1,i2,i3,i4,i5,i6,i7,i8,1)>=-er_th && out_rcr2(i1,i2,i3,i4,i5,i6,i7,i8,1)>=-er_th
                                    i=i+1;                                   
                                    % Final conds with optimization approach
                                    s5=data_out_rcr(i1,i2,i3,i4,i5,i6,i7,i8,1);
                                    tf=data_out_rcr(i1,i2,i3,i4,i5,i6,i7,i8,2);
                                    delta_v=data_out_rcr(i1,i2,i3,i4,i5,i6,i7,i8,3);
                                    v1=[mod(s5, 2*pi) tf delta_v];
                                    
                                    % Final conds with delta_v approach
                                    s52=data_out_rcr2(i1,i2,i3,i4,i5,i6,i7,i8,1);
                                    tf2=data_out_rcr2(i1,i2,i3,i4,i5,i6,i7,i8,2);
                                    delta_v2=wp-we;
                                    v2=[mod(s52, 2*pi) tf2 delta_v2];
                                    
                                    % Error storing
                                    
                                    err_abs(i)=norm(v1-v2);
                                    err_rel(i)=err_abs(i)/norm(v1);
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

figure();
plot(100*err_rel, '-*b');
xlabel('Grid points');
ylabel('Relative distance (%)');
grid on;