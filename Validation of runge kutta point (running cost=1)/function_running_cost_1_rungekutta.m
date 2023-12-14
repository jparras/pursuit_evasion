function [error_norm,error_s5_norm,varargout]=function_running_cost_1_rungekutta(s,param,mode,varargin)
% NOTE: Parameters MUST NOT be normalized!
% Define game parameters and initial conditions
xe0=param.xe0;
ye0=param.ye0;
xp0=param.xp0;
yp0=param.yp0;
up0=param.up0;
vp0=param.vp0;
ue0=param.ue0;
ve0=param.ve0;
wp=param.wp;
Fp=param.Fp;
we=param.we;
Fe=param.Fe;
l=param.l;

ke=Fe/we;
kp=Fp/wp;

s5=s(1);
tf=s(2);

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

% xp(1)=xp0;
% Compute the error norm

dist_final=sqrt((xe(end)-xp(end))^2+(ye(end)-yp(end))^2);
final_angle=atan2(xe(end)-xp(end),ye(end)-yp(end));
if mode==1
    error_norm=dist_final;
    error_s5_norm=abs(mod(final_angle,2*pi)-mod(s(1),2*pi));
    %error_norm=100./(1+exp(-10*(dist_final-l)))+tf;
else
    error_norm=double(dist_final<l);
    error_s5_norm=final_angle;
end

% Additional output
if isempty(varargin)==0
    varargout{1}=tup;
    varargout{2}=tue;
    varargout{3}=tvp;
    varargout{4}=tve;
    varargout{5}=xe;
    varargout{6}=xp;
    varargout{7}=ye;
    varargout{8}=yp;
    varargout{9}=up;
    varargout{10}=ue;
    varargout{11}=vp;
    varargout{12}=ve;
end
        
return;