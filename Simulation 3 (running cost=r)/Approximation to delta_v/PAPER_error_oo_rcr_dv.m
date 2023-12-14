%% ERROR BETWEEN ANALYTICAL EXPRESSION IN GAME WITH L=1 AND OPTIMIZATION METHOD (Appr, Delta_v)
% Juan Parras, GAPS-ETSIT, April 2016
clear all; clc; close all;
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
tf_max=100; %Sufficiently high
Pj=1.11;
Pc=1;
N=100;
A=N*(log2(1+Pc/Pj)+(Pj*(0.0069*D*(1+Pj/Pc)/sqrt(Pj/Pc)+14.4070)/Pc)/(2*D^2*(1+Pj/Pc)^2*log(2)));
B=N*(Pj*log(0.1824*D*(1+Pj/Pc)/sqrt(Pj/Pc)+0.4823)/Pc)/(2*D^2*(1+Pj/Pc)^2*log(2));
dv=wp-we;
denorm_vector=[2*pi tf_max];
global param; % Declared as global, to be used by my_cost_fun
param=struct('xe0',xe0,'ye0',ye0,'xp0',xp0,'yp0',yp0,'up0',up0,'vp0',vp0,'ue0',ue0,'ve0',ve0,'wp',wp,'Fp',Fp,'we',we,'Fe',Fe,'l',l,'ke',ke,'kp',kp,'D',D,'tf_max',tf_max,'A',A,'B',B,'dv',dv);
er_th=0.9;

%% Loop iteration 1
%h=waitbar(0,'Iterating loop ...');
% Out is to store values: 1 for cost on interation, 2 for number of
% iterations used
out=zeros(length(xe0),length(ye0),length(xp0), length(yp0),length(up0), length(vp0), length(ue0), length(ve0),6);
data_out=zeros(length(xe0),length(ye0),length(xp0), length(yp0),length(up0), length(vp0), length(ue0), length(ve0),3);
i=1;
total_grid_length=prod([length(xe0),length(ye0),length(xp0), length(yp0),length(up0), length(vp0), length(ue0), length(ve0)]);
for i1=1:length(xe0)
    param.xe0=xe0(i1);
    for i2=1:length(ye0)
        param.ye0=ye0(i2);
        for i3=1:length(xp0)
            param.xp0=xp0(i3);
            for i4=1:length(yp0)
                param.yp0=yp0(i4);
                for i5=1:length(up0)
                    param.up0=up0(i5);
                    for i6=1:length(vp0)
                        param.vp0=vp0(i6);
                        for i7=1:length(ue0)
                            param.ue0=ue0(i7);
                            for i8=1:length(ve0)
                                param.ve0=ve0(i8);
                                % Optimization in two dimensions: in this
                                % case there is NO analytical solution to
                                % compare with!
                                tic;
                                lp=1;
                                nit=1;
                                while lp
                                    my_cost_fun=@function_running_cost_r_rungekutta_norm;
                                    settings.dim = 2;
                                    settings.type = 'det';
                                    switch nit
                                        case 1
                                            display('1e2 iterations');
                                            x = oo(my_cost_fun ,1e2, settings);
                                            out(i1,i2,i3,i4,i5,i6,i7,i8,2)=1e2;
                                        case 2
                                            display('1e3 iterations');
                                            x = oo(my_cost_fun ,1e3, settings);
                                            out(i1,i2,i3,i4,i5,i6,i7,i8,2)=1e3;
                                        otherwise
                                            display('1e4 iterations');
                                            x = oo(my_cost_fun ,1e4, settings);
                                            out(i1,i2,i3,i4,i5,i6,i7,i8,2)=1e4;
                                    end
                                    x'.*denorm_vector'
                                    function_running_cost_r_rungekutta_norm(x)
                                    
                                    %figure(1);
                                    %plot_traj_param_rcr( x'.*denorm_vector'-[0 0 param.dv/2]', 200)
                                    
                                    out(i1,i2,i3,i4,i5,i6,i7,i8,1)=function_running_cost_r_rungekutta_norm(x); % value function using optimization value
                                    x_oo_denorm=x'.*denorm_vector';
                                    data_out(i1,i2,i3,i4,i5,i6,i7,i8,1)=x_oo_denorm(1);
                                    data_out(i1,i2,i3,i4,i5,i6,i7,i8,2)=x_oo_denorm(2);
                                    if out(i1,i2,i3,i4,i5,i6,i7,i8,1)<=-er_th && nit<3 %No valid solution found with n_iter
                                        lp=1;
                                        nit=nit+1;
                                    else
                                        lp=0; 
                                    end
                                end
                                %waitbar(i/(total_grid_length));
                                save('PAPER_error_oo_data_rcr_dv_appr');
                                toc
                                i=i+1
                            end
                        end
                    end
                end
            end
        end
    end
end
%close(h);
%% Data analytics
pl=1;
if pl
    n_iter=out(:,:,:,:,:,:,:,:,2);
    n_iter=n_iter(out(:,:,:,:,:,:,:,:,1)>=-er_th); %Points correctly solved
    n_iter=n_iter(:);
    figure();
    plot(n_iter); % Number of iterations required to achieve a solution
    figure();
    hist(log10(n_iter));
    

    npor=sum(sum(sum(sum(sum(sum(sum(sum(out(:,:,:,:,:,:,:,:,1)>=-er_th))))))));
    npor1000=sum(n_iter==1000);
    npor10000=sum(n_iter==10000);
    npor100000=sum(n_iter==100000);
    display(['Number of points optimized correctly: ' num2str(npor) ' of ' num2str(total_grid_length) ': ' num2str(100*npor/total_grid_length) ' %']);
    display(['Number of points optimized correctly (1000 iter): ' num2str(npor1000) ' of ' num2str(total_grid_length) ': ' num2str(100*npor1000/total_grid_length) ' %']);
    display(['Number of points optimized correctly (10000 iter): ' num2str(npor10000) ' of ' num2str(total_grid_length) ': ' num2str(100*npor10000/total_grid_length) ' %']);
    display(['Number of points optimized correctly (100000 iter): ' num2str(npor100000) ' of ' num2str(total_grid_length) ': ' num2str(100*npor100000/total_grid_length) ' %']);
end
%% Save
sa=1;
if sa
    save('PAPER_error_oo_data_rcr_dv_appr');
end
display('Simulation finished');