function [flag1, flag2,flag3,flag4,tr_xp,tr_yp,tr_xe,tr_ye,tr_up,tr_vp,tr_ue,tr_ve,xrelay,yrelay,tau,Ct]=main_iso(n_p,xe0,ye0,xp0,yp0,ue0,ve0,up0,vp0,wp,Fp,we,Fe,l,xrelay,yrelay,pc,ij,snr_min,ep,pl)
% varargout:
%function [t_x,t_y,t_xe,t_ye,t_u,t_v,tau_v,t_iter]=main_iso(n_p,n_p_l,kt,sd_ev,sd_pu)
%WARNING: TESTED ON MATLAB R2015b

% Initial conditions
% n_p=50;
% xe0=10;
% ye0=0;
% xp0=0;
% yp0=0;
% up0=10;
% vp0=2;
% ue0=-1;
% ve0=0;
% wp=10;
% Fp=1;
% we=5;
% Fe=5;
% l=1;

%Initial calculations

ke=Fe/we;
kp=Fp/wp;
flag1=0; %Fail in tau calculations
flag2=0; %Fail in initial conditions calculations, numeric method
flag3=0; %Fail in initial conditions calculations, no symbolic method
flag4=0; %Heading angle changed to solve unambiguity

% Tau (final time) is obtained

tau=solve_tau(xe0,ye0,xp0,yp0,ue0,ve0,up0,vp0,l,Fp,kp,Fe,ke);
if isempty(tau)
    flag1=1;
    tr_xp=zeros(n_p,1);
    tr_yp=zeros(n_p,1);
    tr_xe=zeros(n_p,1);
    tr_ye=zeros(n_p,1);
    tr_up=zeros(n_p,1);
    tr_vp=zeros(n_p,1);
    tr_ue=zeros(n_p,1);
    tr_ve=zeros(n_p,1);
    Ct=zeros(n_p,1);
    tau=0;
    return;
end

% Final conditions are obtained from tau

[s1,s2,s3,s4,s5,s6,s7,flag2]=solve_finals(xe0,ye0,xp0,yp0,ue0,ve0,up0,vp0,l,Fp,kp,Fe,ke,tau);
[s1b,s2b,s3b,s4b,cs5,ss5,s6b,s7b,flag3,flag4]=solve_finals_nosym(xe0,ye0,xp0,yp0,ue0,ve0,up0,vp0,l,Fp,kp,Fe,ke,tau);
if flag2==1
    if flag3==1
        tr_xp=zeros(n_p,1);
        tr_yp=zeros(n_p,1);
        tr_xe=zeros(n_p,1);
        tr_ye=zeros(n_p,1);
        tr_up=zeros(n_p,1);
        tr_vp=zeros(n_p,1);
        tr_ue=zeros(n_p,1);
        tr_ve=zeros(n_p,1);
        Ct=zeros(n_p,1);
        return;
    else
        s1=s1b;
        s2=s2b;
        s3=s3b;
        s4=s4b;
        s6=s6b;
        s7=s7b;
        s5=atan2(ss5,cs5); %Uses the four quadrants!
    end
end
% Obtain the trajectories
t0=0;
[tr_xp,tr_yp,tr_xe,tr_ye,tr_up,tr_vp,tr_ue,tr_ve,t]=trajectories(s1,s2,s3,s4,s5,s6,s7,l,Fp,kp,Fe,ke,tau,t0,n_p);

% % Scale relays positions
% 
% xmax=max(max(tr_xe),max(tr_xp));
% xmin=min(min(tr_xe),min(tr_xp));
% ymax=max(max(tr_ye),max(tr_yp));
% ymin=min(min(tr_ye),min(tr_yp));
% % if abs(xmax)>=abs(xmin)
% %     aux=(abs(xmax)+abs(xmin))/2;
% %     xrelay=xrelay*(xmax-aux)+aux;
% % else
% %     aux=(abs(xmax)+abs(xmin))/2;
% %     xrelay=xrelay*(-xmin+aux)+aux;
% % end
% % if abs(ymax)>=abs(ymin)
% %     aux=(abs(ymax)+abs(ymin))/2;
% %     yrelay=yrelay*(ymax-aux)+aux;
% % else
% %     aux=(abs(ymax)+abs(ymin))/2;
% %     yrelay=yrelay*(-ymin+aux)+aux;
% % end
% factor=0.5;
% if xmax>=0
%     xmax=xmax*(1+factor);
% else
%     xmax=xmax/(1+factor);
% end
% if xmin>=0
%     xmin=xmin/(1+factor);
% else
%     xmin=xmin*(1+factor);
% end
% if ymax>=0
%     ymax=ymax*(1+factor);
% else
%     ymax=ymax/(1+factor);
% end
% if ymin>=0
%     ymin=ymin/(1+factor);
% else
%     ymin=ymin*(1+factor);
% end
% 
% xrelay=xrelay*(xmax-xmin)+xmin;
% yrelay=yrelay*(ymax-ymin)+ymin;

%Plot the trajectories

if pl
    
    figure(1);
    plot(tr_xp,tr_yp,'*-r',tr_xe,tr_ye,'o-b',xrelay,yrelay,'*k');grid on;
    title('Trajectories (blue=evader, red=pursuer) and relays (black)');
    xlabel('x');
    ylabel('y');
end
% Check if evader is captured
dist=sqrt((tr_xp-tr_xe).^2+(tr_yp-tr_ye).^2);
if pl
    figure(2)
    subplot(311);
    plot(t(end:-1:1),dist,'*-b',t(end:-1:1),l*ones(1,length(dist)),'*-r');grid on;
    title('Distance between evader and pursuer (blue) and limit of capture (red)');
    xlabel('Time (s)');
    ylabel('Distance (m)');
    subplot(312);
    plot(t(end:-1:1),sqrt(tr_ue.^2+tr_ve.^2),'*-b',t(end:-1:1),we*ones(1,length(t)),'*-k',t(end:-1:1),sqrt(tr_up.^2+tr_vp.^2),'*-r',t(end:-1:1),wp*ones(1,length(t)),'*-g');
    grid on;
    title('Speeds of evader (blue) and its limit (black), and pursuer (red) and its limit (green)');
    xlabel('Time (s)');
    ylabel('Speeds (m/s)');
end
if dist(end)<=l+1e-3 %Threshold for finite precission effects
    sample=find(dist<=1+1e-3,1);
    display(['NO FEEDBACK: capture happened in time = ' num2str(t(end-sample+1),4) ' with tf = ' num2str(t(1),4) ' and final distance = ' num2str(dist(end),4)]);
else
    display(['NO FEEDBACK: capture didn´t happen with final distance = ' num2str(dist(end),4)]);
end

% SNR calculations

SNR=zeros(length(xrelay),n_p);
N0=1e-4;
for i=1:length(xrelay) %For each relay
    for j=1:n_p %For each point in trajectory
        djam=(tr_xp(j)-xrelay(i))^2+(tr_yp(j)-yrelay(i))^2+ep^2;
        dcom=(tr_xe(j)-xrelay(i))^2+(tr_ye(j)-yrelay(i))^2+ep^2;
        %SNR(i,j)=pc/ij*djam/dcom;
        SNR(i,j)=(pc/dcom)/(N0+ij/djam);
    end
end

% Capacity calculations

Ct=zeros(n_p,1); % CApacity in each time point
for i=1:length(Ct)
    aux=SNR(:,i); %SNR in each relay
    aux(aux<=snr_min)=0; % If SNR is below the threshold, is set to zero
    Ct(i)=sum(log2(1+aux));
end

if pl
    subplot(313);
    plot(t(end:-1:1),Ct);
    title(['Evolution of total capacity per bandwidth unit. Final value = ' num2str(Ct(end)) ]);
    xlabel('Time (s)');
    ylabel('SE (bps/Hz)');
end

end

function [tau]=solve_tau(xe,ye,xp,yp,ue,ve,up,vp,l,Fp,kp,Fe,ke)

    syms tau positive;
    eq=(xp-xe-up*(exp(-kp*tau)-1)/kp+ue*(exp(-ke*tau)-1)/ke)^2+(yp-ye-vp*(exp(-kp*tau)-1)/kp+ve*(exp(-ke*tau)-1)/ke)^2-(Fe*(-1+exp(-ke*tau)+ke*tau)/ke^2-l-Fp*(-1+exp(-kp*tau)+kp*tau)/kp^2)^2;
    tau=vpasolve(eq,tau,[0 Inf]);
    tau=double(tau);
end

function [s1,s2,s3,s4,s5,s6,s7,flag]=solve_finals(xe0,ye0,xp0,yp0,ue0,ve0,up0,vp0,l,Fp,kp,Fe,ke,tau)

    syms s1 s2 s3 s4 s5 s6 s7 real
    eq_up=s3*exp(kp*tau)+Fp*sin(s5)*(1-exp(kp*tau))/kp-up0;
    eq_vp=s4*exp(kp*tau)+Fp*cos(s5)*(1-exp(kp*tau))/kp-vp0;
    eq_ue=s6*exp(ke*tau)+Fe*sin(s5)*(1-exp(ke*tau))/ke-ue0;
    eq_ve=s7*exp(ke*tau)+Fe*cos(s5)*(1-exp(ke*tau))/ke-ve0;
    eq_xp=s1+s3*(1-exp(kp*tau))/kp+Fp*sin(s5)*(exp(kp*tau)-1-kp*tau)/kp^2-xp0;
    eq_yp=s2+s4*(1-exp(kp*tau))/kp+Fp*cos(s5)*(exp(kp*tau)-1-kp*tau)/kp^2-yp0;
    eq_xe=s1+l*sin(s5)+s6*(1-exp(ke*tau))/ke+Fe*sin(s5)*(exp(ke*tau)-1-ke*tau)/ke^2-xe0;
    eq_ye=s2+l*cos(s5)+s7*(1-exp(ke*tau))/ke+Fe*cos(s5)*(exp(ke*tau)-1-ke*tau)/ke^2-ye0; 

    [s1,s2,s3,s4,s5,s6,s7]=vpasolve([eq_vp,eq_xp,eq_yp,eq_xe,eq_ye,eq_ue,eq_ve],[s1,s2,s3,s4,s5,s6,s7]);
    s1=double(s1);
    s2=double(s2);
    s3=double(s3);
    s4=double(s4);
    s5=double(s5);
    s6=double(s6);
    s7=double(s7);
    if abs(s3*exp(kp*tau)+Fp*sin(s5)*(1-exp(kp*tau))/kp-up0)<0.01 % chech on up, the equation not used for system solution
        display('Solution found correctly');
        flag=0;
    else
        display(['Systems solution check : ' num2str(s3*exp(kp*tau)+Fp*sin(s5)*(1-exp(kp*tau))/kp-up0)]);
        flag=1;
    end

end

function [s1,s2,s3,s4,cs5,ss5,s6,s7,flag1,flag2]=solve_finals_nosym(xe0,ye0,xp0,yp0,ue0,ve0,up0,vp0,l,Fp,kp,Fe,ke,tau)

    Ap=exp(kp*tau);
    Bp=(1-exp(kp*tau))*Fp/kp;
    Cp=(exp(kp*tau)-1-kp*tau)*Fp/kp^2;
    Dp=(1-exp(kp*tau))/kp;
    Ae=exp(ke*tau);
    Be=(1-exp(ke*tau))*Fe/ke;
    Ce=(exp(ke*tau)-1-ke*tau)*Fe/ke^2+l;
    De=(1-exp(ke*tau))/ke;
    
    s1=(Dp*(-Bp*De*ue0-Ae*Ce*up0+Be*De*up0+Ae*Bp*xe0)+Ap*(Cp*De*ue0-Ae*Cp*xe0+Ae*Ce*xp0-Be*De*xp0))/(-Ap*Be*De+Ae*(Ap*(Ce-Cp)+Bp*Dp));
    s2=(Dp*(-Bp*De*ve0-Ae*Ce*vp0+Be*De*vp0+Ae*Bp*ye0)+Ap*(Cp*De*ve0-Ae*Cp*ye0+Ae*Ce*yp0-Be*De*yp0))/(-Ap*Be*De+Ae*(Ap*(Ce-Cp)+Bp*Dp));
    s3=(Ae*(Ce-Cp)*up0-Be*De*up0+Bp*(De*ue0+Ae*(-xe0+xp0)))/(-Ap*Be*De+Ae*(Ap*(Ce-Cp)+Bp*Dp));
    s4=(Ae*(Ce-Cp)*vp0-Be*De*vp0+Bp*(De*ve0+Ae*(-ye0+yp0)))/(-Ap*Be*De+Ae*(Ap*(Ce-Cp)+Bp*Dp));
    ss5=(-Ae*Dp*up0+Ap*(De*ue0+Ae*(-xe0+xp0)))/(-Ap*Be*De+Ae*(Ap*(Ce-Cp)+Bp*Dp));
    cs5=(-Ae*Dp*vp0+Ap*(De*ve0+Ae*(-ye0+yp0)))/(-Ap*Be*De+Ae*(Ap*(Ce-Cp)+Bp*Dp));
    s6=(Ap*Ce*ue0-Ap*Cp*ue0+Bp*Dp*ue0-Be*Dp*up0-Ap*Be*xe0+Ap*Be*xp0)/(-Ap*Be*De+Ae*(Ap*(Ce-Cp)+Bp*Dp));
    s7=(Ap*Ce*ve0-Ap*Cp*ve0+Bp*Dp*ve0-Be*Dp*vp0-Ap*Be*ye0+Ap*Be*yp0)/(-Ap*Be*De+Ae*(Ap*(Ce-Cp)+Bp*Dp));
    %Sometimes there is ambigüity in s5 calculation: it is s5 or s5+pi
    if abs(ss5^2+cs5^2-1)<0.01 % chech on up, the equation not used for system solution
        display('Solution found correctly');
        flag1=0;
        if abs(s1+s3*(1-exp(kp*tau))/kp+Fp*ss5*(exp(kp*tau)-1-kp*tau)/kp^2-xp0)>=0.01 %Use the correct heading final angle if the one found is not OK
            display('Final heading angle changed');
            cs5=-cs5;
            ss5=-ss5;
            flag2=1;
        end
    else
        display(['Systems solution check : ' num2str(ss5^2+cs5^2-1)]);ss5^2+cs5
        flag1=1;
        flag2=0;
    end

end


function [tr_xp,tr_yp,tr_xe,tr_ye,tr_up,tr_vp,tr_ue,tr_ve,t]=trajectories(s1,s2,s3,s4,s5,s6,s7,l,Fp,kp,Fe,ke,tau,t0,n_p)
    t=tau-linspace(t0,tau,n_p); %Time vector
    tr_up=s3*exp(kp*t)+Fp*sin(s5)*(1-exp(kp*t))/kp;
    tr_vp=s4*exp(kp*t)+Fp*cos(s5)*(1-exp(kp*t))/kp;
    tr_ue=s6*exp(ke*t)+Fe*sin(s5)*(1-exp(ke*t))/ke;
    tr_ve=s7*exp(ke*t)+Fe*cos(s5)*(1-exp(ke*t))/ke;
    tr_xp=s1+s3*(1-exp(kp*t))/kp+Fp*sin(s5)*(exp(kp*t)-1-kp*t)/kp^2;
    tr_yp=s2+s4*(1-exp(kp*t))/kp+Fp*cos(s5)*(exp(kp*t)-1-kp*t)/kp^2;
    tr_xe=s1+l*sin(s5)+s6*(1-exp(ke*t))/ke+Fe*sin(s5)*(exp(ke*t)-1-ke*t)/ke^2;
    tr_ye=s2+l*cos(s5)+s7*(1-exp(ke*t))/ke+Fe*cos(s5)*(exp(ke*t)-1-ke*t)/ke^2; 
end