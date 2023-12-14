%% PURSUIT EVASION OF TWO UAVs WITH IMPULSE CONSTRAINS

% Juan Parras, UPM, January 2016
% Code executed on MatLAB R2015b

clear all;clc;close all;
%% Data inicialization
% Simulation parameters
n_p=100*100; %Number of points per trajectory
pl=1; %To plot data
n_sim=1; %Number of simulations to run
sa=0; % To save data
% Dynamic and capture characteristics
wp=2;
Fp=1;
we=1;
Fe=1;
l=1;
L=100;
% RF characteristics
pc=1;
ij=1.11;
snr_min=1;
nrelays=100;
ep=1;
% Initial position limits
xmax=10; %Max square size, centered at the origin
ymax=10; %Max squared size, centered at the origin

%% Data structure inicialization
fl1=zeros(n_sim,1);
fl2=zeros(n_sim,1);
fl3=zeros(n_sim,1);
fl4=zeros(n_sim,1);
xp=zeros(n_sim,n_p);
yp=zeros(n_sim,n_p);
xe=zeros(n_sim,n_p);
ye=zeros(n_sim,n_p);
up=zeros(n_sim,n_p);
vp=zeros(n_sim,n_p);
ue=zeros(n_sim,n_p);
ve=zeros(n_sim,n_p);
xrel=zeros(n_sim,nrelays);
yrel=zeros(n_sim,nrelays);
ta=zeros(n_sim,1);
C=zeros(n_sim,n_p);
xe0=zeros(n_sim,1);
ye0=zeros(n_sim,1);
xp0=zeros(n_sim,1);
yp0=zeros(n_sim,1);
up0=zeros(n_sim,1);
vp0=zeros(n_sim,1);
ue0=zeros(n_sim,1);
ve0=zeros(n_sim,1);
%% Simulation
for i=1:n_sim
    i
    xe0(i)=(rand(1)-0.5)*xmax;
    ye0(i)=(rand(1)-0.5)*ymax;
    xp0(i)=(rand(1)-0.5)*xmax;
    yp0(i)=(rand(1)-0.5)*ymax;
    up0(i)=rand(1)*wp;
    vp0(i)=rand(1)*sqrt(wp^2-up0(i)^2); %Scaled not to superate wp
    ue0(i)=rand(1)*we;
    ve0(i)=rand(1)*sqrt(we^2-ue0(i)^2); %Scaled not to superate we
    % The relays positions are scaled in the main function!! (Here are
    % normalized)
    xrelay=-L+2*L*rand(nrelays,1); % Values scaled in the function
    yrelay=-L+2*L*rand(nrelays,1); % Values scaled in the function
    xe0=10;
    ye0=0;
    xp0=0;
    yp0=0;
    up0=10;
    vp0=2;
    ue0=-1;
    ve0=0;
    [fl1(i), fl2(i),fl3(i),fl4(i), xp(i,:), yp(i,:), xe(i,:), ye(i,:), up(i,:), vp(i,:), ue(i,:), ve(i,:), xrel(i,:), yrel(i,:), ta(i), C(i,:)]=main_iso(n_p,xe0(i),ye0(i),xp0(i),yp0(i),ue0(i),ve0(i),up0(i),vp0(i),wp,Fp,we,Fe,l,xrelay,yrelay,pc,ij,snr_min,ep,pl);
end
%% Store results
if sa
    save(['Simulation_results_def_paper_short_' num2str(n_sim)]);
end