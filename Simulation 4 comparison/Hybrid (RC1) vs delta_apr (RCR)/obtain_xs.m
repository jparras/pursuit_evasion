function[xs,t]=obtain_xs(s5,tf,delta_v,n_points);
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