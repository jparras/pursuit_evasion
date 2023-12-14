function [] = plot_traj_param_rcr(s,val)
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
A=parameters.A;
B=parameters.B;
l=parameters.l;

ke=Fe/we;
kp=Fp/wp;

s5=s(1);
tf=s(2);
delta_v=s(3);

[t,xs]=ode15s(@(t,x)[...
    x(2);...
    Fp*((eps+2*B*(exp(kp*(tf-t)).*(kp*(tf-t)-1)+1).*(x(5)-x(1))+kp*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(kp*(tf-t))-1).*sin(s5))./sqrt(eps+((eps+2*B*(exp(kp*(tf-t)).*(kp*(tf-t)-1)+1).*(x(5)-x(1))+kp*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(kp*(tf-t))-1).*sin(s5))).^2+((eps+2*B*(exp(kp*(tf-t)).*(kp*(tf-t)-1)+1).*(x(7)-x(3))+kp*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(kp*(tf-t))-1).*cos(s5))).^2))-kp*x(2);...
    x(4);...
    Fp*((eps+2*B*(exp(kp*(tf-t)).*(kp*(tf-t)-1)+1).*(x(7)-x(3))+kp*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(kp*(tf-t))-1).*cos(s5))./sqrt(eps+((eps+2*B*(exp(kp*(tf-t)).*(kp*(tf-t)-1)+1).*(x(5)-x(1))+kp*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(kp*(tf-t))-1).*sin(s5))).^2+((eps+2*B*(exp(kp*(tf-t)).*(kp*(tf-t)-1)+1).*(x(7)-x(3))+kp*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(kp*(tf-t))-1).*cos(s5))).^2))-kp*x(4);...
    x(6);...
    Fe*((eps+2*B*(exp(ke*(tf-t)).*(ke*(tf-t)-1)+1).*(x(5)-x(1))+ke*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(ke*(tf-t))-1).*sin(s5))./sqrt(eps+((eps+2*B*(exp(ke*(tf-t)).*(ke*(tf-t)-1)+1).*(x(5)-x(1))+ke*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(ke*(tf-t))-1).*sin(s5))).^2+((eps+2*B*(exp(ke*(tf-t)).*(ke*(tf-t)-1)+1).*(x(7)-x(3))+ke*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(ke*(tf-t))-1).*cos(s5))).^2))-ke*x(6);...
    x(8);...
    Fe*((eps+2*B*(exp(ke*(tf-t)).*(ke*(tf-t)-1)+1).*(x(7)-x(3))+ke*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(ke*(tf-t))-1).*cos(s5))./sqrt(eps+((eps+2*B*(exp(ke*(tf-t)).*(ke*(tf-t)-1)+1).*(x(5)-x(1))+ke*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(ke*(tf-t))-1).*sin(s5))).^2+((eps+2*B*(exp(ke*(tf-t)).*(ke*(tf-t)-1)+1).*(x(7)-x(3))+ke*((A+B*((x(7)-x(3)).^2+(x(5)-x(1)).^2))./delta_v)*(exp(ke*(tf-t))-1).*cos(s5))).^2))-ke*x(8)],...
    linspace(0,tf,val),[xp0;up0;yp0;vp0;xe0;ue0;ye0;ve0]);
xp=xs(:,1);
up=xs(:,2);
yp=xs(:,3);
vp=xs(:,4);
xe=xs(:,5);
ue=xs(:,6);
ye=xs(:,7);
ve=xs(:,8);

angle=linspace(0,2*pi,50);

%subplot(1,2,1)
% plot(xe,ye,'-b',xp,yp,'--r',xe(1),ye(1),'bo',xp(1),yp(1),'ro',xp(end)+l*cos(angle),yp(end)+l*sin(angle),'k');
% grid on;
% xlabel('X');
% ylabel('Y');
%subplot(1,2,2)
plot(t,sqrt(vp.^2+up.^2),'--r',t,sqrt(ve.^2+ue.^2),'-b',t,wp*ones(1,length(t)),':r',t,we*ones(1,length(t)),':b');
grid on;
xlabel('t');
ylabel('Speed');