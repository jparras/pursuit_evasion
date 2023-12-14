%% ERROR BETWEEN ANALYTICAL EXPRESSION IN GAME WITH L=1 AND OPTIMIZATION METHOD
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
tf_max=[]; % Initialized here, but updated in loop!
global param; % Declared as global, to be used by my_cost_fun
param=struct('xe0',xe0(1),'ye0',ye0(1),'xp0',xp0(1),'yp0',yp0(1),'up0',up0(1),'vp0',vp0(1),'ue0',ue0(1),'ve0',ve0(1),'wp',wp,'Fp',Fp,'we',we,'Fe',Fe,'l',l,'ke',ke,'kp',kp,'D',D,'tf_max',tf_max);
er_th=0.9;

%% Loop iteration 1
%h=waitbar(0,'Iterating loop ...');
% Out is to store values: 1 for abs distance, 2 for rel distance, 3 for
% analytical optimal value, 4 for optimal optimized value, 5 for flag if
% solution found analytically, 6 for number of iterations used
out=zeros(length(xe0),length(ye0),length(xp0), length(yp0),length(up0), length(vp0), length(ue0), length(ve0),6);
data_out=zeros(length(xe0),length(ye0),length(xp0), length(yp0),length(up0), length(vp0), length(ue0), length(ve0),4);
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
                                % Analytical data obtention
                                s_correct_denorm=obtain_analytical_s;
                                if sum(isnan(s_correct_denorm))==0 %Check that there is a valid analytical solution
                                    out(i1,i2,i3,i4,i5,i6,i7,i8,5)=1;
                                end

                                % tf obtention
                                if out(i1,i2,i3,i4,i5,i6,i7,i8,5)==1
                                    param.tf_max=ceil(s_correct_denorm(2)/10)*10;
                                else
                                    param.tf_max=4*D; %Pretty high value
                                end
                                denorm_vector=[2*pi param.tf_max];
                                s_correct_norm=s_correct_denorm./denorm_vector;
                                % Optimization in two dimensions
                                opt_x =s_correct_norm';
                                lp=1;
                                nit=1;
                                while lp
                                    my_cost_fun = @function_running_cost_1_rungekutta_norm;
                                    settings.dim = 2;
                                    settings.type = 'det';
                                    switch nit
                                        case 1
                                            display('1e3 iterations');
                                            x = oo(my_cost_fun ,1e3, settings);
                                            out(i1,i2,i3,i4,i5,i6,i7,i8,6)=1e3;
                                        case 2
                                            display('1e4 iterations');
                                            x = oo(my_cost_fun ,1e4, settings);
                                            out(i1,i2,i3,i4,i5,i6,i7,i8,6)=1e4;
                                        otherwise
                                            display('1e5 iterations');
                                            x = oo(my_cost_fun ,1e5, settings);
                                            out(i1,i2,i3,i4,i5,i6,i7,i8,6)=1e5;
                                    end
                                    [opt_x.*denorm_vector' x'.*denorm_vector']
                                    if out(i1,i2,i3,i4,i5,i6,i7,i8,5)==1
                                        [function_running_cost_1_rungekutta_norm(opt_x') function_running_cost_1_rungekutta_norm(x)]
                                    else
                                        function_running_cost_1_rungekutta_norm(x)
                                    end
                                    out(i1,i2,i3,i4,i5,i6,i7,i8,1)=norm(opt_x.*denorm_vector'-x'.*denorm_vector');
                                    out(i1,i2,i3,i4,i5,i6,i7,i8,2)=out(i1,i2,i3,i4,i5,i6,i7,i8,1)/norm(opt_x.*denorm_vector'); %Normalized with respect analytical value norm vector
                                    if out(i1,i2,i3,i4,i5,i6,i7,i8,5)==1
                                        out(i1,i2,i3,i4,i5,i6,i7,i8,3)=function_running_cost_1_rungekutta_norm(opt_x'); % value function using analytical value
                                    else
                                        out(i1,i2,i3,i4,i5,i6,i7,i8,3)=100; % Value for incorrect cost! It's a flag!
                                    end
                                    out(i1,i2,i3,i4,i5,i6,i7,i8,4)=function_running_cost_1_rungekutta_norm(x); % value function using optimization value
                                    x_oo_denorm=x'.*denorm_vector';
                                    data_out(i1,i2,i3,i4,i5,i6,i7,i8,1)=s_correct_denorm(1);
                                    data_out(i1,i2,i3,i4,i5,i6,i7,i8,2)=s_correct_denorm(2);
                                    data_out(i1,i2,i3,i4,i5,i6,i7,i8,3)=x_oo_denorm(1);
                                    data_out(i1,i2,i3,i4,i5,i6,i7,i8,4)=x_oo_denorm(2);
                                    if out(i1,i2,i3,i4,i5,i6,i7,i8,4)<=-er_th && nit<3 %No valid solution found with n_iter
                                        lp=1;
                                        nit=nit+1;
                                    else
                                        lp=0; 
                                    end
                                end
                                %waitbar(i/(total_grid_length));
                                save('PAPER_error_oo_analytical_data_rc1');
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
pl=0;
if pl
    n_iter=out(:,:,:,:,:,:,:,:,6);
    n_iter=n_iter(out(:,:,:,:,:,:,:,:,4)>=-er_th); %Points correctly solved
    n_iter=n_iter(:);
    figure();
    plot(n_iter); % Number of iterations required to achieve a solution
    figure();
    hist(log10(n_iter));
    rel_dist=out(:,:,:,:,:,:,:,:,2);
    rel_dist=rel_dist(out(:,:,:,:,:,:,:,:,4)>=-er_th); %Points correctly solved
    rel_dist=rel_dist(:);
    figure();
    plot(rel_dist*100,'*-b');
    grid on;
    xlabel('Grid point');
    ylabel('Relative distance (%)');


    npor=sum(sum(sum(sum(sum(sum(sum(sum(out(:,:,:,:,:,:,:,:,4)>=-er_th))))))));
    npor1000=sum(n_iter==1000);
    npor10000=sum(n_iter==10000);
    npor100000=sum(n_iter==100000);
    npar=sum(sum(sum(sum(sum(sum(sum(sum(out(:,:,:,:,:,:,:,:,5)>=er_th))))))));
    display(['Number of points optimized correctly: ' num2str(npor) ' of ' num2str(total_grid_length) ': ' num2str(100*npor/total_grid_length) ' %']);
    display(['Number of points optimized correctly (1000 iter): ' num2str(npor1000) ' of ' num2str(total_grid_length) ': ' num2str(100*npor1000/total_grid_length) ' %']);
    display(['Number of points optimized correctly (10000 iter): ' num2str(npor10000) ' of ' num2str(total_grid_length) ': ' num2str(100*npor10000/total_grid_length) ' %']);
    display(['Number of points optimized correctly (100000 iter): ' num2str(npor100000) ' of ' num2str(total_grid_length) ': ' num2str(100*npor100000/total_grid_length) ' %']);
    display(['Number of points solved analytically correctly: ' num2str(npar) ' of ' num2str(total_grid_length) ': ' num2str(100*npar/total_grid_length) ' %']);
end
%% Save
sa=1;
if sa
    save('PAPER_error_oo_analytical_data_rc1');
end
display('simulation finished');