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
global param; % Declared as global, to be used by my_cost_fun
param=struct('xe0',xe0,'ye0',ye0,'xp0',xp0,'yp0',yp0,'up0',up0,'vp0',vp0,'ue0',ue0,'ve0',ve0,'wp',wp,'Fp',Fp,'we',we,'Fe',Fe,'l',l,'ke',ke,'kp',kp,'D',D,'tf_max',tf_max);

%% Analytical data obtention
s_correct_denorm=obtain_analytical_s;
param.tf_max=ceil(s_correct_denorm(2)/10)*10;
denorm_vector=[2*pi param.tf_max];
s_correct_norm=s_correct_denorm./denorm_vector;

%% Optimization
opt_x =s_correct_norm';

my_cost_fun = @function_running_cost_1_rungekutta_norm;
settings.dim = 2;
settings.type = 'det';
x = oo(my_cost_fun ,1e3, settings);

%[opt_x, x']
[opt_x.*denorm_vector' x'.*denorm_vector']
[function_running_cost_1_rungekutta_norm(opt_x') function_running_cost_1_rungekutta_norm(x)]
% norm(opt_x - x')
% [my_cost_fun(opt_x'), my_cost_fun(x)]

%% Plotting of value function
pl=0;
if pl
    s5=linspace(0,2*pi,40)/(2*pi);
    tf=linspace(0,param.tf_max,40)/param.tf_max;
    f=0;
    for i=1:length(s5)
        for j=1:length(tf)
            f(i,j)=function_running_cost_1_rungekutta_norm([s5(i) tf(j)]);
        end
    end
    mesh(s5*2*pi,tf*param.tf_max,f');
    xlabel('s5');
    ylabel('tf');
end
