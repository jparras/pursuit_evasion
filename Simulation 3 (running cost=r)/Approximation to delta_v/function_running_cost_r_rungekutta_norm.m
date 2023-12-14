function [error_norm]=function_running_cost_r_rungekutta_norm(s)

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
tf_max=parameters.tf_max;
ke=parameters.ke;
kp=parameters.kp;
delta_v=parameters.dv;
A=parameters.A;
B=parameters.B;

% Denormalization
denorm_vector=[2*pi tf_max];
s=s.*denorm_vector;
 
s5=s(1);
tf=s(2);


% Compute initial values of trajectories with provided final conditions
% Use of Runge-Kutta method to solve the differential equation (2nd order!)

if tf==0 %If tf==0, ode45 throws an error
    tf=tf+2*eps;
end
% Notation: x(1)=xp, x(3)=yp, x(5)=xe, x(7)=ye, x(2)=up, x(4)=vp, x(6)=ue,
% x(8)=ve
w='on';
if strcmp(w,'off')
    warning('off','MATLAB:illConditionedMatrix');
    warning('off','MATLAB:ode15s:IntegrationTolNotMet');
else
    warning('on','MATLAB:illConditionedMatrix');
    warning('on','MATLAB:ode15s:IntegrationTolNotMet');
end
[t,xs]=ode15s(@(t,x)[...
    x(2);...
    Fp*((eps+2*B*(exp(kp*(tf-t)).*(kp*(tf-t)-1)+1).*(x(5)-x(1))+kp*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(kp*(tf-t))-1).*sin(s5))./sqrt(eps+((eps+2*B*(exp(kp*(tf-t)).*(kp*(tf-t)-1)+1).*(x(5)-x(1))+kp*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(kp*(tf-t))-1).*sin(s5))).^2+((eps+2*B*(exp(kp*(tf-t)).*(kp*(tf-t)-1)+1).*(x(7)-x(3))+kp*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(kp*(tf-t))-1).*cos(s5))).^2))-kp*x(2);...
    x(4);...
    Fp*((eps+2*B*(exp(kp*(tf-t)).*(kp*(tf-t)-1)+1).*(x(7)-x(3))+kp*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(kp*(tf-t))-1).*cos(s5))./sqrt(eps+((eps+2*B*(exp(kp*(tf-t)).*(kp*(tf-t)-1)+1).*(x(5)-x(1))+kp*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(kp*(tf-t))-1).*sin(s5))).^2+((eps+2*B*(exp(kp*(tf-t)).*(kp*(tf-t)-1)+1).*(x(7)-x(3))+kp*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(kp*(tf-t))-1).*cos(s5))).^2))-kp*x(4);...
    x(6);...
    Fe*((eps+2*B*(exp(ke*(tf-t)).*(ke*(tf-t)-1)+1).*(x(5)-x(1))+ke*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(ke*(tf-t))-1).*sin(s5))./sqrt(eps+((eps+2*B*(exp(ke*(tf-t)).*(ke*(tf-t)-1)+1).*(x(5)-x(1))+ke*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(ke*(tf-t))-1).*sin(s5))).^2+((eps+2*B*(exp(ke*(tf-t)).*(ke*(tf-t)-1)+1).*(x(7)-x(3))+ke*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(ke*(tf-t))-1).*cos(s5))).^2))-ke*x(6);...
    x(8);...
    Fe*((eps+2*B*(exp(ke*(tf-t)).*(ke*(tf-t)-1)+1).*(x(7)-x(3))+ke*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(ke*(tf-t))-1).*cos(s5))./sqrt(eps+((eps+2*B*(exp(ke*(tf-t)).*(ke*(tf-t)-1)+1).*(x(5)-x(1))+ke*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(ke*(tf-t))-1).*sin(s5))).^2+((eps+2*B*(exp(ke*(tf-t)).*(ke*(tf-t)-1)+1).*(x(7)-x(3))+ke*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(ke*(tf-t))-1).*cos(s5))).^2))-ke*x(8)],...
    [eps tf],[xp0;up0;yp0;vp0;xe0;ue0;ye0;ve0]);
xp=xs(:,1);
up=xs(:,2);
yp=xs(:,3);
vp=xs(:,4);
xe=xs(:,5);
ue=xs(:,6);
ye=xs(:,7);
ve=xs(:,8);

% Compute the error norm
dist_final=sqrt((xe(end)-xp(end))^2+(ye(end)-yp(end))^2);
final_angle=atan2(xe(end)-xp(end),ye(end)-yp(end));
v_f_p=sqrt(up(end)^2+vp(end)^2);
v_f_e=sqrt(ue(end)^2+ve(end)^2);
delta_v_final=v_f_p-v_f_e;

error_norm= - (1./(1+exp(-500*(dist_final-l)))+abs(mod(final_angle,2*pi)-mod(s5,2*pi))+abs(delta_v_final-delta_v));
if isnan(error_norm)
    error_norm=-10;
end
%error_norm= -dist_final;

return;