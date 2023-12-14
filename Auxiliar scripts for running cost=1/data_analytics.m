%% PURSUIT EVASION WITH IMPULSE CONSTRAINS DATA ANALYSIS TOOL
% Juan Parras, GAPS - UPM, January 2016
clear all;clc; close all;
%% Load data
uiopen('matlab');

%% Data analytics
% Flag errors calculations
display(['Errors due to tau calculations (%)= ' num2str(sum(fl1)*100/n_sim)]);
display(['Errors due to system solver (%) = ' num2str(sum(fl2)*100/n_sim)]);
display(['Errors due to system solver not "fixed" (%) = ' num2str(sum(fl2.*fl3)*100/n_sim)]);
display(['Heading angle changes (%) = ' num2str(sum(fl4)*100/n_sim)]);
% Error calculations
er_tot=0;
th=0.01; %Threshold margin
for i=1:n_sim
    % Check that initial conditions are satisfied in all cases
    er_xp=abs(xp(i,1)-xp0(i));
    er_yp=abs(yp(i,1)-yp0(i));
    er_xe=abs(xe(i,1)-xe0(i));
    er_ye=abs(ye(i,1)-ye0(i));
    er_up=abs(up(i,1)-up0(i));
    er_vp=abs(vp(i,1)-vp0(i));
    er_ue=abs(ue(i,1)-ue0(i));
    er_ve=abs(ve(i,1)-ve0(i));
    if er_xp>th || er_yp>th || er_xe>th || er_ye>th || er_up>th || er_vp>th || er_ue>th || er_ve>th
        er_tot=er_tot+1;
    end
end
display(['TOTAL ERRORS (%) = ' num2str(er_tot*100/n_sim)]);
    
% Trajectories plot
opt=input('¿Desea ver todas las trayectorias (1) o sólo las correctas (2)?');
if opt==1
    v=fl2.*fl3;
else
    v=fl2;
end

ang=0:0.01:2*pi; 
xcirc=l*cos(ang);
ycirc=l*sin(ang);

figure(1);

for i=1:n_sim
    if fl1(i)==0&&v(i)==0
        % Data preparation
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
        % Plot
        subplot(3,4,[1:3 5:7 9:11]);
        plot(tr_xp,tr_yp,'*-r',tr_xe,tr_ye,'o-b',xrelay,yrelay,'*k',tr_xe(end)+xcirc,tr_ye(end)+ycirc,'g',xp0(i)+xcirc,yp0(i)+ycirc,'m',[xp0(i),xp0(i)+up0(i)*2*l],[yp0(i),yp0(i)+vp0(i)*2*l],'m',xe0(i)+xcirc,ye0(i)+ycirc,'m',[xe0(i),xe0(i)+ue0(i)*2*l],[ye0(i),ye0(i)+ve0(i)*2*l],'m');grid on;
        title({'Trajectories (blue=evader, red=pursuer) and relays (black).Capture circle is in green. Initial consitions in magent.',['Flag in symbolic system solution = ' num2str(fl2(i)) ';Flag in no symbolic solution = ' num2str(fl3(i))]});
        xlabel('x');
        ylabel('y');
        subplot(3,4,4);
        plot(t(end:-1:1),dist,'*-b',t(end:-1:1),l*ones(1,length(t)),'*-r');grid on;
        title({'Distance between evader and pursuer (blue)','and limit of capture (red)'});
        xlabel('Time (s)');
        ylabel('Distance (m)');
        subplot(3,4,8);
        plot(t(end:-1:1),sqrt(tr_ue.^2+tr_ve.^2),'*-b',t(end:-1:1),we*ones(1,length(t)),'*-k',t(end:-1:1),sqrt(tr_up.^2+tr_vp.^2),'*-r',t(end:-1:1),wp*ones(1,length(t)),'*-g');
        grid on;
        title({'Speeds of evader (blue) and its limit (black)','and pursuer (red) and its limit (green)'});
        xlabel('Time (s)');
        ylabel('Speeds (m/s)');
        subplot(3,4,12);
        plot(t(end:-1:1),Ct);
        title({'Evolution of total capacity per bandwidth unit',['Final value = ' num2str(Ct(end)) ]});
        xlabel('Time (s)');
        ylabel('SE (bps/Hz)');
        pause();
    end
end