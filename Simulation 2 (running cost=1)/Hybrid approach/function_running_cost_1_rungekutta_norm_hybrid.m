function [error_norm]=function_running_cost_1_rungekutta_norm_hybrid(s)

parameters=getGlobalx;
xe0=parameters.xe0;
ye0=parameters.ye0;
xp0=parameters.xp0;
yp0=parameters.yp0;
up0=parameters.up0;
vp0=parameters.vp0;
ue0=parameters.ue0;
ve0=parameters.ve0;
wp=parameters.wp;
Fp=parameters.Fp;
we=parameters.we;
Fe=parameters.Fe;
l=parameters.l;
ke=parameters.ke;
kp=parameters.kp;
tf=parameters.tf;

% Denormalization
s=s*2*pi;
 
s5=s;

% Compute initial values of trajectories with provided final conditions
% Use of Runge-Kutta method to solve the differential equation (2nd order!)

if tf==0 %If tf==0, ode45 throws an error
    tf=tf+eps;
end

% [tup,up]=ode45(@(t,up) Fp*sin(s5)-kp*up,[0 tf],up0);
% [tvp,vp]=ode45(@(t,vp) Fp*cos(s5)-kp*vp,[0 tf],vp0);
% [tue,ue]=ode45(@(t,ue) Fe*sin(s5)-ke*ue,[0 tf],ue0);
% [tve,ve]=ode45(@(t,ve) Fe*cos(s5)-ke*ve,[0 tf],ve0);
options = odeset('RelTol',1e-10);
[tup,xs1]=ode45(@(t,x1)[x1(2);Fp*sin(s5)-kp*x1(2)],[0 tf],[xp0;up0],options);
xp=xs1(:,1);
up=xs1(:,2);
[tvp,xs2]=ode45(@(t,x2)[x2(2);Fp*cos(s5)-kp*x2(2)],[0 tf],[yp0;vp0],options);
yp=xs2(:,1);
vp=xs2(:,2);
[tue,xs3]=ode45(@(t,x3)[x3(2);Fe*sin(s5)-ke*x3(2)],[0 tf],[xe0;ue0],options);
xe=xs3(:,1);
ue=xs3(:,2);
[tve,xs4]=ode45(@(t,x4)[x4(2);Fe*cos(s5)-ke*x4(2)],[0 tf],[ye0;ve0],options);
ye=xs4(:,1);
ve=xs4(:,2);

% Compute the error norm
dist_final=sqrt((xe(end)-xp(end))^2+(ye(end)-yp(end))^2);
final_angle=atan2(xe(end)-xp(end),ye(end)-yp(end));

error_norm= - (1./(1+exp(-500*(dist_final-l)))+abs(mod(final_angle,2*pi)-mod(s5,2*pi)));
%error_norm= -dist_final;

return;