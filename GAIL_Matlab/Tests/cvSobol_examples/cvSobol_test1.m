%  test cubSobol and cvSobol on simple func 

% init and set number of iteration
clc; clearvars;
iter=5; 
t=0; n=0;   t1=0;t2=0;t3=0;t4=0;    n1=0; n2=0;n3=0;n4=0;

%% define functions and parameters
% simple eg
%f=@(x) x + 0.1*x.^2; g=@(x) x-0.5; d=1; 
%f=@(x)sum(x+0.1*x.^2,2); g=@(x) sum(x-0.5,2); d=4; 
%f=@(x) sin( (pi/2).*x ); g=@(x) -x.^2+2.*x-2/3; d=1; 

% multi dim example 1: prod
%f=@(x) x(:,1).^2 .*x(:,2).^2 .* x(:,3).^2; g=@(x) x(:,1) .*x(:,2) .* x(:,3)-0.5^3;d=3;

% multi dim example 2: sum
%f=@(x) sin( (pi/2).*x(:,1) )+ sin( (pi/2).*x(:,2) ) + sin( (pi/2).*x(:,3) );
%g=@(x) (-x(:,1).^2+2.*x(:,1)-2/3) + (-x(:,2).^2+2.*x(:,2)-2/3)+ (-x(:,3).^2+2.*x(:,3)-2/3);d=3;

% multi cv example 1: dim=1

%f=@(x) x.^3; g={@(x) x-1/2, @(x) x.^2- 1/3}; d=1; 

% multi cv example 2: dim > 1 

f=@(x) x(:,1)+ x(:,2).^2 +0.01*x(:,3).^3; g={@(x) x(:,1)-1/2, @(x) x(:,2).^2- 1/3}; d=3; 
% comparison between multi cv and single cv
%f=@(x) x(:,1)+ x(:,2).^2 +0.01*x(:,3).^3; g=@(x) (x(:,1)-1/2)+ (x(:,2).^2- 1/3); d=3; 
 
r_lag=0;
abstol=1e-6;
fprintf('\n abstol: %d \n', abstol);
%% cubSobol
% run it a couple times first to stablize it 
for i=1:5
	cvSobol_a0(f,d,abstol,0,'uniform');
end

% begin testing
for i=1:iter
	[q,out]=cvSobol_a0(f,d,abstol,0,'uniform');
	t=t+out.time;
    n=n+out.n;
end
fprintf('\n Results of CubSobol_g: \n');
fprintf('q=%.10f \n',q);
fprintf('avg time of cubSobol: %s \n', num2str(t/iter) ); 
fprintf('avg n of cubSobol: %s \n', num2str(n/iter) );


%% cvSobol_a1

% run it a couple times first to make it stable 
for i=1:5
	cvSobol_a1(f,g,d,abstol,0,'uniform');
end

% begin testing
for i=1:iter
	[q1,out1]=cvSobol_a1(f,g,d,abstol,0,'uniform');
	t1=t1+out1.time;
    n1=n1+out1.n;
end
fprintf('\n Results of cvSobol_a1(L2 Reg): \n');
fprintf('q1=%.10f  \n',q1);
fprintf('q1-q=%d  \n',q1-q);
fprintf('avg time of cvSobol_a1: %s \n', num2str(t1/iter) ); 
fprintf('avg n of cvSobol_a1: %s \n', num2str(n1/iter) );

%{
%% cvSobol_a2

% run it a couple times first to elimiate start up error
for i=1:5
	cvSobol_a1(f,g,d,abstol,0,'uniform',r_lag, 'L2', 'T');
end

% begin testing
for i=1:iter
	[q2,out2]=cvSobol_a1(f,g,d,abstol,0,'uniform',r_lag, 'L2', 'T');
	t2=t2+out2.time;
    n2=n2+out2.n;
end
fprintf('\n Results of cvSobol_a2(L2 Reg with updates): \n');
fprintf('q2=%d  \n',q2);
fprintf('avg time of cvSobol_a2: %s \n', num2str(t2/iter) ); 
fprintf('avg n of cvSobol_a2: %s \n', num2str(n2/iter) );


%% cvSobol_a3

% run it a couple times first to elimiate start up error
for i=1:5
	cvSobol_a1(f,g,d,abstol,0,'uniform', r_lag, 'L1');
end

% begin testing
for i=1:iter
	[q3,out3]=cvSobol_a1(f,g,d,abstol,0,'uniform', r_lag, 'L1');
	t3=t3+out3.time;
    n3=n3+out3.n;
end
fprintf('\n Results of cvSobol_a3(L1 Reg): \n');
fprintf('q3=%d  \n',q3);
fprintf('avg time of cvSobol_a3: %s \n', num2str(t3/iter) ); 
fprintf('avg n of cvSobol_a3: %s \n', num2str(n3/iter) );


%% cvSobol_a4

% run it a couple times first to elimiate start up error
for i=1:5
	cvSobol_a4(f,g,d,abstol,0,'uniform', r_lag);
end

% begin testing
for i=1:iter
	[q4,out4]=cvSobol_a4(f,g,d,abstol,0,'uniform', r_lag);
	t4=t4+out4.time;
    n4=n4+out4.n;
end
fprintf('\n Results of cvSobol_a4(L1 Reg with updates): \n');
fprintf('q4=%d  \n',q4);
fprintf('avg time of cvSobol_a4: %s \n', num2str(t4/iter) ); 
fprintf('avg n of cvSobol_a4: %s \n', num2str(n4/iter) );
%}
