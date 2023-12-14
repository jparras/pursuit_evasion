function [s_correct_denorm]=obtain_analytical_s(param);

%% Input data
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
D=param.D; %Maximum abs value of the coordinates

%% tf calculation
syms tau positive;
eq=(xp0-xe0-up0*(exp(-kp*tau)-1)/kp+ue0*(exp(-ke*tau)-1)/ke)^2+(yp0-ye0-vp0*(exp(-kp*tau)-1)/kp+ve0*(exp(-ke*tau)-1)/ke)^2-(Fe*(-1+exp(-ke*tau)+ke*tau)/ke^2-l-Fp*(-1+exp(-kp*tau)+kp*tau)/kp^2)^2;
tau=vpasolve(eq,tau, [0 Inf]); %Numerical solution is faster: no analytical solution for expression anyway
tau=double(tau);

%% s5 calculation
flag=1;
if isempty(tau)==0
    flag=0;
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
        if abs(s1+s3*(1-exp(kp*tau))/kp+Fp*ss5*(exp(kp*tau)-1-kp*tau)/kp^2-xp0)>=0.01 %Use the correct heading final angle if the one found is not OK
            display('Final heading angle changed');
            cs5=-cs5;
            ss5=-ss5;
        end
    else
        display(['Systems solution check : ' num2str(ss5^2+cs5^2-1)]);
        flag=1;
    end
    s5=atan2(ss5,cs5);
    s5=mod(s5,2*pi);
end
%% Return data
if flag==0
    if tau>0
        s_correct_denorm=[s5 tau];
    else
        s_correct_denorm=[NaN NaN]; %Solution found is not right: tf must be positive!
    end
else
    s_correct_denorm=[NaN NaN]; %Solution found is not right
end
return;
