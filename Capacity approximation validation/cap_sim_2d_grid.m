%% Capacity simulation in 2D

%Juan Parras, GAPS-UPM, February 2016
clear all;close all; clc;
%% Parameters setup
pl=1; %To plot data
sa=0; %To save data
n_sim=10; %Number of simulations to run
L=100;
delta=0.5; %Space step
xc=10;
yc=0;
xj=-L:delta:L;
yj=xj;
N=100; %Number of relays
ep=1; 

%% Calculations
C=zeros(length(xj),length(yj),n_sim);
for ind=1:n_sim
    ind
    xr=-L+2*L*rand(N,1);
    yr=-L+2*L*rand(N,1);
    for i=1:N %For each relay
        dcom(i)=(xc-xr(i))^2+(yc-yr(i))^2+ep^2;
    end
    for j=1:length(xj) %For all jamm xj
        for k=1:length(yj) %For all jamm yj
            SNR=zeros(1,N);
            for l=1:N %For each relay
                djam=(xj(j)-xr(l))^2+(yj(k)-yr(l))^2+ep^2;
                SNR(l)=djam/dcom(l);
            end
            C_aux(j,k)=sum(log2(1+SNR));
        end
    end
    C(:,:,ind)=C_aux;
end
%% Theoretical results (using approximation)
k=1; %Pj=Pc considered in above section
for i=1:length(xj) %For all jamm xj
    for j=1:length(yj) %For all jamm yj
        A=k*((yc-yj(j)).^2+(xc-xj(i)).^2).*asinh(L*(1+k)./(sqrt(k*((yc-yj(j)).^2+(xc-xj(i)).^2))));
        if isnan(A)
            A=0;% Taking the limit when xc=xj,yc=yj
        end
        C_teor(i,j)=N*(log2(1+1/k)+A/(2*L^2*(1+k)^2*log(2)));
    end
end
%% Minimum search

Cm=mean(C,3);
[~,i2]=min(min(Cm,[],1));
[~,i1]=min(min(Cm,[],2));
d_emp=sqrt((xj(i1)-xc)^2+(yj(i2)-yc)^2);
[~,i2t]=min(min(C_teor,[],1));
[~,i1t]=min(min(C_teor,[],2));
d_teor=sqrt((xj(i1t)-xc)^2+(yj(i2t)-yc)^2);

if pl
%     figure();
%     contour(xj,yj,Cm);
%     hold on;
%     plot(xc,yc,'b*',xj(i1),yj(i2),'ro');
%     title('Theoretical contours of Capacity obtained');
%     xlabel('xj');
%     ylabel('yj');
%     hold off;
    figure();
    mesh(xj,yj,Cm);
    hold on;
    mesh(xj,yj,C_teor);
    title('Capacity functions obtained: theoretical vs empirical');
    xlabel('xj');
    ylabel('yj');
    hold off;
    figure();
    imagesc(xj,yj,Cm-C_teor); colorbar;
    title('Difference of capacity functions obtained - VALIDATES THE APPROXIMATION!');
    xlabel('xj');
    ylabel('yj');
end
%% Error graph
% In this section, an error graph is created. With communicator in (xc,yc),
% the jammer is placed over the main diagonal (xj=yj) and the theoretical
% capacity value is computed. The difference between that value and the
% empirical will be the absolute error, whereas relative will be divided by
% the empirical result
er_rel=zeros(1,length(xj));
er_abs=zeros(1,length(xj));
for i=1:length(xj)
%     A=k*((yc-yj(i)).^2+(xc-xj(i)).^2).*asinh(L*(1+k)./(sqrt(k*((yc-yj(i)).^2+(xc-xj(i)).^2))));
%     if isnan(A)
%         A=0;% Taking the limit when xc=xj,yc=yj
%     end
%     C_t(i)=N*(log2(1+1/k)+A/(2*L^2*(1+k)^2*log(2)));
    er_abs(i)=abs(C_teor(i,i)-Cm(i,i));
    er_rel(i)=er_abs(i)/Cm(i,i);
end
figure();
plot(xj,er_abs,'*-');
xlabel('x_j=y_j');
ylabel('Absolute error (bps/Hz)');
grid on;
figure();
plot(xj,er_rel*100,'*-');
xlabel('x_j=y_j');
ylabel('Relative error (%)');
grid on;
%% Theoretical results, without approximation of the function
for i=1:length(xj) %For all jamm xj
    for j=1:length(yj) %For all jamm yj
        A1=((1+2*k+k^2)*L+(1+k)*yj(j)+(k+k^2)*yc);
        A2=sqrt((1+2*k+k^2)*L^2+((2*k^2+2*k)*yc+(2*k+2)*yj(j))*L+(1+k)*yj(j).^2+(k+k^2)*yc^2+k*xj(i)^2-2*k*xc*xj(i)+k*xc^2);
        A3=sqrt(k^2+2*k+1)*(k*yj(j)^2-2*k*yc*yj(j)+k*yc^2+k*xj(i)^2-2*k*xc*xj(i)+k*xc^2)*asinh((sqrt((4*k+8*k^2+4*k^3)*yj(j)^2+(-8*k-16*k^2-8*k^3)*yc*yj(j)+(4*k+8*k^2+4*k^3)*yc^2+(4*k+8*k^2+4*k^3)*xj(i)^2+(-8*k-16*k^2-8*k^3)*xc*xj(i)+(4*k+8*k^2+4*k^3)*xc^2)*(k*yc+yj(j)+(k+1)*L))/((2*k+2*k^2)*yj(j)^2+(-4*k-4*k^2)*yc*yj(j)+(2*k+2*k^2)*yc^2+(2*k+2*k^2)*xj(i)^2+(-4*k-4*k^2)*xc*xj(i)+(2*k+2*k^2)*xc^2));
        A4=((1+2*k+k^2)*L+(-1-k)*yj(j)+(-k-k^2)*yc);
        A5=sqrt((1+2*k+k^2)*L^2+((-2*k^2-2*k)*yc+(-2*k-2)*yj(j))*L+(1+k)*yj(j)^2+(k+k^2)*yc^2+k*xj(i)^2-2*k*xc*xj(i)+k*xc^2);
        A6=sqrt(k^2+2*k+1)*(k*yj(j)^2-2*k*yc*yj(j)+k*yc^2+k*xj(i)^2-2*k*xc*xj(i)+k*xc^2)*asinh((sqrt((4*k+8*k^2+4*k^3)*yj(j)^2+(-8*k-16*k^2-8*k^3)*yc*yj(j)+(4*k+8*k^2+4*k^3)*yc^2+(4*k+8*k^2+4*k^3)*xj(i)^2+(-8*k-16*k^2-8*k^3)*xc*xj(i)+(4*k+8*k^2+4*k^3)*xc^2)*(-k*yc-yj(j)+(k+1)*L))/((2*k+2*k^2)*yj(j)^2+(-4*k-4*k^2)*yc*yj(j)+(2*k+2*k^2)*yc^2+(2*k+2*k^2)*xj(i)^2+(-4*k-4*k^2)*xc*xj(i)+(2*k+2*k^2)*xc^2));
        if isnan(A3)
            A3=0;% Taking the limit when xc=xj,yc=yj
        end
        if isnan(A6)
            A6=0;% Taking the limit when xc=xj,yc=yj
        end
        A=(A1*A2+A3)/(2*k^2+4*k+2)+(A4*A5+A6)/(2*k^2+4*k+2);
        C_teor_napr(i,j)=N*pi/(2*L^2*(1+k)*log(2))*A+((2*log((k+1)/k)-pi)*N)/(2*log(2));
    end
end
%% Minimum search

Cm=mean(C,3);
[~,i2]=min(min(Cm,[],1));
[~,i1]=min(min(Cm,[],2));
d_emp=sqrt((xj(i1)-xc)^2+(yj(i2)-yc)^2);
[~,i2t]=min(min(C_teor_napr,[],1));
[~,i1t]=min(min(C_teor_napr,[],2));
d_teor=sqrt((xj(i1t)-xc)^2+(yj(i2t)-yc)^2);

if pl
%     figure();
%     contour(xj,yj,Cm);
%     hold on;
%     plot(xc,yc,'b*',xj(i1),yj(i2),'ro');
%     title('Theoretical contours of Capacity obtained');
%     xlabel('xj');
%     ylabel('yj');
%     hold off;
    figure();
    mesh(xj,yj,Cm);
    hold on;
    mesh(xj,yj,C_teor_napr);
    title('Capacity functions obtained: theoretical vs empirical');
    xlabel('xj');
    ylabel('yj');
    hold off;
    figure();
    imagesc(xj,yj,Cm-C_teor_napr); colorbar;
    title('Difference of capacity functions obtained - VALIDATES THE APPROXIMATION!');
    xlabel('xj');
    ylabel('yj');
end
%% Error graph

er_rel_napr=zeros(1,length(xj));
er_abs_napr=zeros(1,length(xj));
for i=1:length(xj)
    er_abs_napr(i)=abs(C_teor_napr(i,i)-Cm(i,i));
    er_rel_napr(i)=er_abs_napr(i)/Cm(i,i);
end
figure();
plot(xj,er_abs_napr,'*-');
xlabel('x_j=y_j');
ylabel('Absolute error (bps/Hz)');
grid on;
figure();
plot(xj,er_rel_napr*100,'*-');
xlabel('x_j=y_j');
ylabel('Relative error (%)');
grid on;

%% Comparative error graph

figure();
s=plot(xj,er_abs,'-b',xj,er_abs_napr,'--r');
set(s,'LineWidth',2);
xlabel('x_j=y_j');
ylabel('Absolute error (bps/Hz)');
grid on;
figure();
s=plot(xj,er_rel*100,'-b',xj,er_rel_napr*100,'--r');
set(s,'LineWidth',2);
xlabel('x_j=y_j');
ylabel('Relative error (%)');
grid on;
%% Store results
if sa
    save(['2dcap_sim_results_' num2str(n_sim)]);
end


                
                