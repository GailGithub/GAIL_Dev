clc;clearvars;
iter=5;
t=0; n=0;   t1=0;t2=0;t3=0;t4=0;    n1=0; n2=0;n3=0;n4=0;
randtype='SOBOL'; %IID, SOBOL, or LATTICE sampling
path_param.disctype='BB'; %timestep or BB time differencing
path_param.d=8; %number of trading periods
path_param.T=8/52; %time to expiry
path_param.S0=100; %initial stock price
path_param.sig=0.5; %volatility
path_param.r=0.015; %interest rate
path_param.meanshift=false; %mean shift importance sampling
path_param.drift=-3; %drift in mean
path_param.anti=false; %antithetic variates
pay_param.K=110; %strike price
pay_param.barrier=250; %barrier price

% param for control variates
path_param1.disctype='BB'; %timestep or BB time differencing
path_param1.d=8; %number of trading periods
path_param1.T=8/52; %time to expiry
path_param1.S0=100; %initial stock price
path_param1.sig=0.5; %volatility
path_param1.r=0.015; %interest rate
path_param1.meanshift=false; %mean shift importance sampling
path_param1.drift=-3; %drift in mean
path_param1.anti=false; %antithetic variates
pay_param1.K=110; %strike price
pay_param1.barrier=250; %barrier price

% param for control variates 2
path_param2.disctype='BB'; %timestep or BB time differencing
path_param2.d=8; %number of trading periods
path_param2.T=8/52; %time to expiry
path_param2.S0=100; %initial stock price
path_param2.sig=0.5; %volatility
path_param2.r=0.015; %interest rate
path_param2.meanshift=false; %mean shift importance sampling
path_param2.drift=-3; %drift in mean
path_param2.anti=false; %antithetic variates
pay_param2.K=110; %strike price
pay_param2.barrier=250; %barrier price

r_lag=8;
abstol=1e-4; %absolute error tolerance
fprintf('\n abstol=%d \n', abstol);

%% define target func 
%pay_param.paytype=['amerput']; name='amerput';
pay_param.paytype=['ameanput']; name='ameanput';
pay_param.control={}; %no control variates
f=@(x) payoffx(x, pay_param, path_param);

%hyperbox=[-inf(1,path_param.d);inf(1,path_param.d)];

%% define cv func 
%pay_param1.paytype=['europut']; name='europut';
pay_param1.paytype=['gmeanput']; name1='gmeanput';
pay_param1.control={}; %no control variates
pay_param2.paytype=['gmeancall']; name2='gmeancall';
pay_param2.control={}; %no control variates


% exact price of options chosen as cv 
%put1=exactoptionprice(path_param1, pay_param1, pay_param1.paytype); put1=put1(1);
%g=@(x) payoffx(x, pay_param1,path_param1)-put1;

% multi cv example
put=exactoptionprice(path_param1, pay_param1, pay_param1.paytype); put=put(1);
call=exactoptionprice(path_param2, pay_param2, pay_param2.paytype); call=call(1);
g={@(x) payoffx(x, pay_param1,path_param1)-put, @(x) payoffx(x, pay_param2, path_param2)-call};


% run it a couple times first to elimiate start up error
for i=1:5
	cvSobol_a0(f,path_param.d, abstol, 0,'normal');
end

% begin testing
for i=1:iter
	[q,out]=cvSobol_a0(f,path_param.d, abstol, 0,'normal');
	t=t+out.time;
    n=n+out.n;
end
fprintf('\n Results of CubSobol_g: \n');
fprintf('q=%8.5f  \n',q);
fprintf('avg time of cubSobol: %s \n', num2str(t/iter) ); 
fprintf('avg n of cubSobol: %s \n', num2str(n/iter) );

%% cvSobol_a1

% run it a couple times first to elimiate start up error
for i=1:5
	cvSobol_a1(f,g,path_param.d, abstol,0,'normal', r_lag);
end

% begin testing
for i=1:iter
	[q1,out1]=cvSobol_a1(f,g,path_param.d,abstol, 0,'normal', r_lag);
	t1=t1+out1.time;
    n1=n1+out1.n;
end
fprintf('\n Results of cvSobol_a1(L2 reg): \n');
fprintf('q1=%8.5f  \n',q1);
fprintf('q1-q=%.7f  \n',q1-q);
fprintf('avg time of cvSobol_a1: %s \n', num2str(t1/iter) ); 
fprintf('avg n of cvSobol_a1: %s \n', num2str(n1/iter) );


%% cvSobol_a2

% run it a couple times first to elimiate start up error
for i=1:5
	cvSobol_a1(f,g,path_param.d, abstol,0,'normal', r_lag, 'L2', 'T');
end

% begin testing
for i=1:iter
	[q2,out2]=cvSobol_a1(f,g,path_param.d, abstol,0,'normal', r_lag, 'L2', 'T');
	t2=t2+out2.time;
    n2=n2+out2.n;
end
fprintf('\n Results of cvSobol_a2(L2 with updates): \n');
fprintf('q2=%8.5f  \n',q2);
fprintf('q2-q1=%f, q2-q=%f,  \n',q2-q1, q2-q);
fprintf('avg time of cvSobol_a2: %s \n', num2str(t2/iter) ); 
fprintf('avg n of cvSobol_a2: %s \n', num2str(n2/iter) );

%{
%% cvSobol_a3

% run it a couple times first to elimiate start up error
for i=1:5
	cvSobol_a3(f,g,path_param.d, abstol,0,'normal', r_lag);
end

% begin testing
for i=1:iter
	[q3,out3]=cvSobol_a3(f,g,path_param.d, abstol,0,'normal', r_lag);
	t3=t3+out3.time;
    n3=n3+out3.n;
end
fprintf('\n Results of cvSobol_a3(L1 reg): \n');
fprintf('q3=%8.5f  \n',q3);
fprintf('avg time of cvSobol_a3: %s \n', num2str(t3/iter) ); 
fprintf('avg n of cvSobol_a3: %s \n', num2str(n3/iter) );


%% cvSobol_a4

% run it a couple times first to elimiate start up error
for i=1:5
	cvSobol_a4(f,g,path_param.d, abstol,0,'normal', r_lag);
end

% begin testing
for i=1:iter
	[q4,out4]=cvSobol_a4(f,g,path_param.d, abstol,0,'normal', r_lag);
	t4=t4+out4.time;
    n4=n4+out4.n;
end
fprintf('\n Results of cvSobol_a4(L1 with updates): \n');
fprintf('q4=%8.5f  \n',q4);
fprintf('avg time of cvSobol_a4: %s \n', num2str(t4/iter) ); 
fprintf('avg n of cvSobol_a4: %s \n', num2str(n4/iter) );
%}
