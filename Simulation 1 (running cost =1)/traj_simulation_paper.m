i=1;  
% Simulation parameters
n_p=50; %Number of points per trajectory
pl=0; %To plot data
% Dynamic and capture characteristics
wp=4;
Fp=0.05;
we=1;
Fe=0.03;
l=0.1;
L=100;
% RF characteristics
pc=1;
ij=1.01;
snr_min=1;
nrelays=100;
ep=1;
% Initial position limits
xmax=10; %Max square size, centered at the origin
ymax=10; %Max squared size, centered at the origin

xe0(i)=-20;
ye0(i)=20;
xp0(i)=20;
yp0(i)=-20;
up0(i)=-2;
vp0(i)=-2; %Scaled not to superate wp
ue0(i)=0; %-1,0
ve0(i)=-1; %Scaled not to superate we
% The relays positions are scaled in the main function!! (Here are
% normalized)
xrelay=-L+2*L*rand(nrelays,1); % Values scaled in the function
yrelay=-L+2*L*rand(nrelays,1); % Values scaled in the function

[fl1(i), fl2(i),fl3(i),fl4(i), xp(i,:), yp(i,:), xe(i,:), ye(i,:), up(i,:), vp(i,:), ue(i,:), ve(i,:), xrel(i,:), yrel(i,:), ta(i), C(i,:)]=main_iso(n_p,xe0(i),ye0(i),xp0(i),yp0(i),ue0(i),ve0(i),up0(i),vp0(i),wp,Fp,we,Fe,l,xrelay,yrelay,pc,ij,snr_min,ep,pl);

%% Data preparation
tr_xp=xp(i,:);
tr_yp=yp(i,:);
tr_xe=xe(i,:);
tr_ye=ye(i,:);
xrelay=xrel(i,:);
yrelay=yrel(i,:);
tr_up=up(i,:);
tr_vp=vp(i,:);
tr_ue=ue(i,:);
tr_ve=ve(i,:);
Ct=C(i,:);
t=linspace(ta(i),0,n_p);
dist=sqrt((tr_xp-tr_xe).^2+(tr_yp-tr_ye).^2);
%% Plot
figure();
s=plot(tr_xp,tr_yp,'--r',tr_xe,tr_ye,'-b',xrelay,yrelay,'ok');grid on;
set(s,'LineWidth',2);
xlabel('x');
ylabel('y');
%axis([-15 0 -10 20]);
figure();
s=plot(t(end:-1:1),sqrt(tr_ue.^2+tr_ve.^2),'-b',t(end:-1:1),sqrt(tr_up.^2+tr_vp.^2),'--r');
set(s,'LineWidth',2);
grid on;
xlabel('Time (s)');
ylabel('Speeds (m/s)');
figure();
s=plot(t(end:-1:1),C);
set(s,'LineWidth',2);
grid on;
xlabel('Time (s)');
ylabel('Total system capacity (bps/Hz)');