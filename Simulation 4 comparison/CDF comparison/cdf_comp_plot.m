%% CDF plot comparison
%Juan Parras, GAPS-UPM, March 2016
clear all; clc; close all;
%% Load data
load('CDF_analytical');
load('CDF_delta');
%% Plot
th=2;
df=10; %Decimation factor
id1=find(x_e_analytical>=th,1);
id2=find(x_p_analytical>=th,1);
id3=find(x_e_delta>=th,1);
id4=find(x_p_delta>=th,1);
x_e_analytical=x_e_analytical(1:10*df:id1);
y_e_analytical=y_e_analytical(1:10*df:id1);
x_p_analytical=x_p_analytical(1:10*df:id2);
y_p_analytical=y_p_analytical(1:10*df:id2);
x_e_delta=x_e_delta(1:df:id3);
y_e_delta=y_e_delta(1:df:id3);
x_p_delta=x_p_delta(1:df:id4);
y_p_delta=y_p_delta(1:df:id4);
figure();
subplot(1,2,1);
plot(x_e_analytical,y_e_analytical,'b-',x_e_delta,y_e_delta,'r--');
grid on;
xlim([0 th]);
xlabel('\zeta');
ylabel('CDF');
subplot(1,2,2);
plot(x_p_analytical,y_p_analytical,'b-',x_p_delta,y_p_delta,'r--');
grid on;
xlim([0 th]);
xlabel('\zeta');
ylabel('CDF');