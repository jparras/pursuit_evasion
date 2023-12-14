%% ARCSINH APPROXIMATION AS A LINE
% Juan Parras, GAPS-UPM, March 2016
clear all; clc; close all;
%% Parameter introduction

D=100;
n_points=50;
m=linspace(0,5,n_points);
b=linspace(-5,45,n_points);
C=linspace(0,500,n_points);
%% Loop and calculations
h=waitbar(0, 'Wait for the loop...');
out=zeros(length(m),length(b), length(C));
for i=1:length(m)
    for j=1:length(b)
        for k=1:length(C)
            out(i,j,k)=integral(@(r) (m(i).*r+b(j)-r.*asinh(C(k)./sqrt(r))).^2,0,D);
        end
    end
    waitbar(i/length(m));
end
close(h);
%save('Asinh_appr_data')
%% Minima search
result=zeros(length(C),3);
for i=1:length(C) %For each value of C
    result(i,1)=C(i); %Value of C
    aux(:,:)=out(:,:,i);
    [val,ind]=min(aux(:)); % Minimum search
    [A,B]=ind2sub(size(aux),ind);
    result(i,2)=m(A); % Value of slope
    result(i,3)=b(B); %Value of origin ordinate
end

%% Paramenters adjustement
% Using polyfit, degree 1, exp(result(:,2)) gives the following slope adjustement
% Note: C=x is the variable
% Line adjustement
p1=polyfit(result(:,1),exp(result(:,2)),1)
p2=polyfit(result(:,1),result(:,3),1)
figure();
plot(result(:,1),result(:,2),'*b-.',result(:,1),log(p1(1)*result(:,1)+p1(2)),'o--r');
xlabel('K');
ylabel('m');
legend('Grid','Adjusted');
grid on;

figure();
plot(result(:,1),result(:,3),'*b-.',result(:,1),p2(1)*result(:,1)+p2(2),'o--r')
xlabel('K');
ylabel('b');
legend('Grid','Adjusted');
grid on;

%% Error evolution as a function of C
error=zeros(length(C),1);
for i=1:length(C)
    abs_error=sqrt(integral(@(r) (result(i,2).*r+result(i,3)-r.*asinh(result(i,1)./sqrt(r))).^2,0,D)); %Abs error is the sqrt: the area of error
    normalization=integral(@(r) r.*asinh(result(i,1)./sqrt(r)),0,D); % We normalize with the area under the curve
    error(i)=abs_error/normalization;
end
error2=zeros(length(C),1);
for i=1:length(C)
    abs_error=sqrt(integral(@(r) (log(p1(1).*result(i,1)+p1(2)).*r+p2(1).*result(i,1)+p2(2)-r.*asinh(result(i,1)./sqrt(r))).^2,0,D)); %Abs error is the sqrt: the area of error
    normalization=integral(@(r) r.*asinh(result(i,1)./sqrt(r)),0,D); % We normalize with the area under the curve
    error2(i)=abs_error/normalization;
end
plot(C,error*100,'*b-.',C,error2*100,'o--r');
xlabel('K');
ylabel('Relative error (%)');
legend('Grid','Adjusted');
grid on;

%% Graphical checkout
r=linspace(0,D,1e3);
figure();
for i=1:length(C)
    f1=result(i,2)*r+result(i,3); % Line optimal
    f2=log(p1(1)*result(i,1)+p1(2)).*r+p2(1)*result(i,1)+p2(2); %Line adjusted
    f3=r.*asinh(result(i,1)./sqrt(r)); %Asinh
    plot(r,f1,'b-.',r,f2,'--r',r,f3,'-k');
    xlabel('r');
    ylabel('f_1(r)');
    legend('Grid','Adjusted','Original')
    grid on;
    pause();
end    
