%define test function
f1 = @(x) x.^2;
f2 = @(x) exp(-1./(x-0.5).^2);
% f3 = @(x,a) exp(-10000*(x-a).^2);
% f4 = @(x,a) 1/4*a*exp(-2*x).*(a-2*exp(x).*(-1 + a*cos(x) - a*sin(x))...
%     +exp(2*x).*(a + 2*cos(x) - 2* sin(x) - a *sin(2*x)));
f5 = @(x) sin(pi*x);
n = 5;
nrep = 200;
time = zeros(n,2,nrep);
npoints = zeros(n,2,nrep);
a = -2;
b = 2;
abstol = 1e-8;
randa = rand(nrep,1)*(b-a)+a;
warning('off','MATLAB:funappx_g:peaky')
for i = 1:nrep;
    tic;
    [~, out_param] = funappx_g(f1,a,b,abstol,9,100);
    t=toc;
    time(1,1,i) = t;
    npoints(1,1,i) = out_param.npoints;
    tic;
    [~, out_param] = funappxglobal_g(f1,a,b,abstol,9,100);
    t=toc;
    time(1,2,i) = t;
    npoints(1,2,i) = out_param.npoints;
    tic;
    [~, out_param] = funappx_g(f2,a,b,abstol,9,100);
    t=toc;
    time(2,1,i) =  t;
    npoints(2,1,i) = out_param.npoints;
    tic;
    [~, out_param] = funappxglobal_g(f2,a,b,abstol,9,100);
    t=toc;
    time(2,2,i) =  t;
    npoints(2,2,i) = out_param.npoints;
    f3 = @(x) exp(-1000*(x-randa(i)).^2);
    tic;
    [~, out_param] = funappx_g(f3,a,b,abstol,9,100);
    t=toc;
    time(3,1,i) =   t;
    npoints(3,1,i) = out_param.npoints;
    tic;
    [~, out_param] = funappxglobal_g(f3,a,b,abstol,9,100);
    t=toc;
    time(3,2,i) =   t;
    npoints(3,2,i) = out_param.npoints;
    f4 = @(x) 1/4*randa(i)*exp(-2*x).*(randa(i)-2*exp(x).*(-1 +...
        randa(i)*cos(x) - randa(i)*sin(x))+exp(2*x).*(randa(i) + 2*cos(x)...
        - 2* sin(x) - randa(i)*sin(2*x)));
    tic;
    [~, out_param] = funappx_g(f4,a,b,abstol,9,100);
    t=toc;
    time(4,1,i) = t;
    npoints(4,1,i) = out_param.npoints;
    tic;
    [~, out_param] = funappxglobal_g(f4,a,b,abstol,9,100);
    t=toc;
    time(4,2,i) =   t;
    npoints(4,2,i) = out_param.npoints;
    tic;
    [~, out_param] = funappx_g(f5,a,b,abstol,9,100);
    t=toc;
    time(5,1,i) =  t;
    npoints(5,1,i) = out_param.npoints;
    tic;
    [~, out_param] = funappxglobal_g(f5,a,b,abstol,9,100);
    t=toc;
    time(5,2,i) =   t;
    npoints(5,2,i) = out_param.npoints;
end;


warning('on','MATLAB:funappx_g:peaky')
display(' ')
display('   Test        # Points          Time Used')
display(' Function   Local    Global     Local    Global')
for i=1:n
    display(sprintf(['%8.0f %8.0f  %7.0f %10.6f  %10.6f'],...
        [i npoints(i,1) npoints(i,2) time(i,1)/nrep time(i,2)/nrep])) 
end

for i=1:nrep, for j=1:5,  timeratio(i,j) = time(j,1, i)/time(j,2, i); end, end
for i=1:nrep, for j=1:5,  npointsratio(i,j) = npoints(j,1, i)/npoints(j,2, i); end, end

subplot(1,2,1); plot(sort(timeratio(:)))
subplot(1,2,2); plot(sort(npointsratio(:)))


