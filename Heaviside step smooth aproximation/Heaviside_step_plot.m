%% HEAVISIDE SMOOTH APPROXIMATION SCRIPT
%Juan Parras, GAPS-UPM, March 2016

%% Define parameters

t=linspace(-1/2,1/2,1e3);
k_2=[1 5 10 50 100 500];
color=cell(1,6);
color{1}='-b';
color{2}='--r';
color{3}=':g';
color{4}='-.k';
color{5}='-y';
color{6}='--m';
color2=cell(1,6);
color2{1}='-ob';
color2{2}='--*r';
color2{3}=':^g';
color2{4}='-.>k';
color2{5}='-<y';
color2{6}='--xm';
color3=cell(1,6);
color3{1}='ob';
color3{2}='*r';
color3{3}='^g';
color3{4}='>k';
color3{5}='<y';
color3{6}='xm';

%% Plot
figure();
hold on;
for i=1:length(k_2)
    f=1./(1+exp(-k_2(i)*t));
    plot(t(1:100:end),f(1:100:end),color2{i});
end
hold off;
for i=1:length(k_2)
    f=1./(1+exp(-k_2(i)*t));
    plot(t(1:100:end),f(1:100:end),color3{i});
    hold on;
end
for i=1:length(k_2)
    f=1./(1+exp(-k_2(i)*t));
    plot(t,f,color{i});
end
grid;
xlabel('X');
ylabel('Heaviside step function approximations')
legend('1', '5', '10', '50', '100', '500');