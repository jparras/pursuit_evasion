function [] = plot_traj_param_cap(s,val)

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

ke=Fe/we;
kp=Fp/wp;

s5=s(1);
tf=s(2);

options = odeset('RelTol',1e-10);
[tup,xs1]=ode45(@(t,x1)[x1(2);Fp*sin(s5)-kp*x1(2)],linspace(0,tf,val),[xp0;up0],options);
xp=xs1(:,1);
up=xs1(:,2);
[tvp,xs2]=ode45(@(t,x2)[x2(2);Fp*cos(s5)-kp*x2(2)],linspace(0,tf,val),[yp0;vp0],options);
yp=xs2(:,1);
vp=xs2(:,2);
[tue,xs3]=ode45(@(t,x3)[x3(2);Fe*sin(s5)-ke*x3(2)],linspace(0,tf,val),[xe0;ue0],options);
xe=xs3(:,1);
ue=xs3(:,2);
[tve,xs4]=ode45(@(t,x4)[x4(2);Fe*cos(s5)-ke*x4(2)],linspace(0,tf,val),[ye0;ve0],options);
ye=xs4(:,1);
ve=xs4(:,2);
angle=linspace(0,2*pi,50);
%subplot(1,2,1);
plot(xe,ye,'*-b',xp,yp,'*-r',xe(1),ye(1),'bo',xp(1),yp(1),'ro',xp(end)+l*cos(angle),yp(end)+l*sin(angle),'k');
%subplot(1,2,2);
%plot(tup, sqrt(up.^2+vp.^2),'-r',tup, sqrt(ue.^2+ve.^2),'-b');