%Example of Asian Arithmetic Mean Option

%% Garbage collection and initialization
format compact %remove blank lines from output
clear all %clear all variables
close all %close all figures
set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLineMarkerSize',40) %latex axis labels

%% Initialization
path_param.S0=100; %initial stock price
path_param.T=1; %time to expiry
path_param.d=26; %number of time steps
path_param.r=0.015; %interest rate
path_param.sig=0.4; %volatility
pay_param.K=110; %strike price
%pay_param.paytype='ameancall'; name='Arithmetic Mean Asian Call'; 
pay_param.paytype='gmeancall'; name='Geometric Mean Asian Call'; 
    %payoff type
abstol=0.1; %tolerance

%% Computation
pay=@(n) payoff(n,pay_param,path_param);
    %function to generate payoffs
[price,out_param]=meanMC_g(pay,'abstol',abstol);
    %compute price

%% Output
disp('For a stock modeled by the geometric Brownian motion')
disp(['initial stock price = $' num2str(path_param.S0)])
disp(['       time horizon = ' num2str(path_param.T) ' years'])
disp(['      interest rate = ' num2str(100*path_param.r) '%'])
disp(['         volatility = ' num2str(path_param.sig)])
disp(['For the ' name ' Option'])
disp(['       # time steps = ' num2str(path_param.d)])
disp(['  approximate price = $' num2str(price)])
disp(['          tolerance = $' num2str(out_param.abstol)])
disp(['        uncertainty = ' num2str(100*out_param.alpha) '%'])
disp(['      total samples = ' num2str(out_param.n)])
disp(['       time elapsed = ' num2str(out_param.time) ' seconds' char(10)])


