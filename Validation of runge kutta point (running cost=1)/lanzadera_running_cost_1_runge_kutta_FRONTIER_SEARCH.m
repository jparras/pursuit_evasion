%% SCRIPT FOR FINDING FRONTIER USING RUNGE-KUTTA WITH RUNNING COST=1
% Juan Parras, GAPS-UPM, March 2016
close all; clc; clear all;
%% Parameter inicialization
n_p=29;
mode=1; % Error is final distance between players
n_it=10;

xe0=10;
ye0=10;
xp0=-10;
yp0=-10;
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
%tf_max=4*D/(wp-we); %High enough

param=struct('xe0',xe0,'ye0',ye0,'xp0',xp0,'yp0',yp0,'up0',up0,'vp0',vp0,'ue0',ue0,'ve0',ve0,'wp',wp,'Fp',Fp,'we',we,'Fe',Fe,'l',l,'ke',ke,'kp',kp,'D',D);
%% Loop iterations
error_vector=zeros(n_it,1);
final_time_1=zeros(n_it,1);
final_time_2=zeros(n_it,1);
final_time_3=zeros(n_it,1);


for k=1:n_it
    display(['Iteration ' num2str(k) ' of ' num2str(n_it)]);
    if n_it~=1
        initial_cond=1;
        while initial_cond
            % Random initial conditions generation, bounded by max allowed
            param.xe0=xe0*rand(1);
            param.ye0=ye0*rand(1);
            param.ue0=ue0*rand(1);
            param.ve0=ve0*rand(1);
            param.xp0=xp0*rand(1);
            param.yp0=yp0*rand(1);
            param.up0=up0*rand(1);
            param.vp0=vp0*rand(1);
            param.wp=wp*rand(1);
            param.Fp=Fp*rand(1);
            param.we=we*rand(1);
            param.Fe=Fe*rand(1);
            param.ke=param.Fe/param.we;
            param.kp=param.Fp/param.wp;
            % Analytical data obtention
            s_correct_denorm=obtain_analytical_s(param);
            if sum(isnan(s_correct_denorm))==0 %Check that there is a valid solution, otherwise, change initial conditions
                initial_cond=0;
            end
        end
    else
        s_correct_denorm=obtain_analytical_s(param);
    end
    param
    s_correct_denorm
    % denorm_vector=[2*pi tf_max];
    % s_correct_norm=s_correct_denorm./denorm_vector;

    % Loop iteration and error storing
    n_p_aux=ceil(n_p/2);
    s5a=linspace(s_correct_denorm(1)*0.8,s_correct_denorm(1),n_p_aux);
    s5b=linspace(s_correct_denorm(1),s_correct_denorm(1)/0.8,n_p_aux);
    s5=[s5a s5b(2:end)];
    tfa=linspace(s_correct_denorm(2)*0.8,s_correct_denorm(2),n_p_aux);
    tfb=linspace(s_correct_denorm(2),s_correct_denorm(2)/0.8,n_p_aux);
    tf=[tfa tfb(2:end)];    
    error=zeros(n_p,n_p);
    error_s5=zeros(n_p,n_p);
    h=waitbar(0,'Computing loop...');
    for i=1:n_p
    for j=1:n_p
        [error(i,j) error_s5(i,j)]=function_running_cost_1_rungekutta([s5(i) tf(j)],param,mode); %Error is final distance between players!
    end
    waitbar(i/n_p);
    end
    close(h);


    % figure();
    % surf(s5,tf,error');
    % title('Error (Final distance) obtained');
    % xlabel('s_5');
    % ylabel('t_f');
    
    %display(['Minimum distance obtained: ' num2str(min(min(error)))]);
    
%     mat_aux=zeros(n_p,n_p);
%     mat_aux(ceil(n_p/2),ceil(n_p/2))=0.5;
%     
%     figure();
%     imagesc(s5,tf,double(error<l)+mat_aux);colorbar;
%     title('Region frontier obtained and analytical solution (dot)');
%     xlabel('s_5');
%     ylabel('t_f');
%   
    if error(ceil(n_p/2),ceil(n_p/2))>(l+1e-3) %Analytical and symbolic solution do not agree: a threshold of 1e-3 is used to avoid numerical issues
        error_vector(k)=1;
        save(['Error_values_in_iteration_' num2str(k)]);
        display('ERROR FOUND!!')
    end
    % Relative error calculation
    % In this game, payoff is final time. We will store final time using
    % three different criteria: analytical solution, minimum norm solution
    % and minimun norm of angle solution, in order to compare
    % the error between them
    
    %Analytical solution
    final_time_1(k)=s_correct_denorm(2);
    % Minimum distance norm solution
    val_min=min(min(error));
    index=find(error==val_min);
    [A B]=ind2sub(size(error),index);
    final_time_2(k)=tf(B);
    % Minimum s5 norm
    val_min=min(min(error_s5(error<l)));
    index=find(error_s5==val_min);
    [A B]=ind2sub(size(error_s5),index);
    final_time_3(k)=tf(B);
    %pause();
end
figure();
rel_error=abs(final_time_1-final_time_2)./final_time_1;
plot(1:n_it,rel_error*100,'b*-');
title('Relative error in each iteration, using distance norm');
xlabel('Iteration');
ylabel('Rel error (%)');
grid on;

figure();
rel_error=abs(final_time_1-final_time_3)./final_time_1;
plot(1:n_it,rel_error*100,'b*-');
title('Relative error in each iteration, using final angle norm');
xlabel('Iteration');
ylabel('Rel error (%)');
grid on;

flags=(error<l);
if max(max(flag))==0
    error('Not enough points for grid');
end
