clear variables
clc
close all

% %% definition of the function to be optimized
% sin1 = @(value) (sin(13 * value) * sin(27 * value) / 2.0 + 0.5); % used in ICML 2013 paper
% guirland =  @(x) 4*x*(1-x)*(0.75+0.25*(1-sqrt(abs(sin(60*x))))); % used in ICML 2013 paper
% 
% difficult =  @(x) 1-sqrt(x) + (-x*x +sqrt(x) )*(sin(1/(x*x*x))+1)/2;
% 
% %% select which function from above we want to optimize
% myfun = sin1; fmax = 0.975599143811574975870826165191829204559326171875;
% % myfun = guirland; fmax = 0.997772313413222;
% % myfun = difficult; 
% 
% %% definitions of the auxiliary functions based on myfun
% myfun_minus = @(value) -myfun(value);
% myfun_noise = @(value) myfun(value) + (rand-0.5)/10;
% myfun_noise_minus = @(value) -myfun_noise(value);

% calling stosoo with 1000 evaluations of myfun_noise
%% Parameters initialization
xe0=10;
ye0=0;
xp0=0;
yp0=0;
up0=10;
vp0=2;
ue0=-1;
ve0=0;
wp=2;
Fp=1;
we=1;
Fe=1;
l=1;
ke=Fe/we;
kp=Fp/wp;
D=100; %Maximum abs value of the coordinates
tf_max=2*D/(wp-we); %High enough
Pj=1.11;
Pc=1;
D=100;
N=100;
A=N*(log2(1+Pc/Pj)+(Pj*(0.0069*D*(1+Pj/Pc)/sqrt(Pj/Pc)+14.4070)/Pc)/(2*D^2*(1+Pj/Pc)^2*log(2)));
B=N*(Pj*log(0.1824*D*(1+Pj/Pc)/sqrt(Pj/Pc)+0.4823)/Pc)/(2*D^2*(1+Pj/Pc)^2*log(2));
dv=wp-we;
denorm_vector=[2*pi tf_max dv];
global param; % Declared as global, to be used by my_cost_fun
param=struct('xe0',xe0,'ye0',ye0,'xp0',xp0,'yp0',yp0,'up0',up0,'vp0',vp0,'ue0',ue0,'ve0',ve0,'wp',wp,'Fp',Fp,'we',we,'Fe',Fe,'l',l,'ke',ke,'kp',kp,'D',D,'tf_max',tf_max,'A',A,'B',B,'dv',dv);

%% Optimization

my_cost_fun = @function_running_cost_r_rungekutta_norm;
settings.dim = 3;
settings.type = 'det';
x = oo(my_cost_fun ,1e3, settings);

%[opt_x, x']
x'.*denorm_vector'
function_running_cost_r_rungekutta_norm(x)
% norm(opt_x - x')
% [my_cost_fun(opt_x'), my_cost_fun(x)]
