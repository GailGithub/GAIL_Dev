%comparison between funappx_g and funappxlocal_g
nrep = 100;
c = rand(100,1)*4;
abstol = 1e-7;
n = 4;
npoints = zeros(n,2,nrep);
time = zeros(n,2,nrep);
warning('off','MATLAB:funappx_g:peaky')
warning('off','MATLAB:funappxglobal_g:peaky')
warning('off','MATLAB:funappxglobal_g:exceedbudget')
for i = 1:nrep;
    a = -c(i) -1;
    b = c(i)+1;
    f1 = @(x) c(i)*x.^2;
    f2 = @(x) sin(c(i)*pi*x);
    f3 = @(x) exp(-1000*(x-c(i)).^2);
    f4 = @(x) 1/4*c(i)*exp(-2*x).*(c(i)-2*exp(x).*(-1 +...
        c(i)*cos(x) - c(i)*sin(x))+exp(2*x).*(c(i) + 2*cos(x)...
        - 2* sin(x) - c(i)*sin(2*x)));
    tic;
    [~, out_param] = funappx_g(f1,a,b,abstol,100,1000);
    t=toc;
    time(1,1,i) = t;
    npoints(1,1,i) = out_param.npoints;
    tic;
    [~, out_param] = funappxglobal_g(f1,a,b,abstol,100,1000);
    t=toc;
    time(1,2,i) = t;
    npoints(1,2,i) = out_param.npoints;
    tic;
    [~, out_param] = funappx_g(f2,a,b,abstol,100,1000);
    t=toc;
    time(2,1,i) =  t;
    npoints(2,1,i) = out_param.npoints;
    tic;
    [~, out_param] = funappxglobal_g(f2,a,b,abstol,100,1000);
    t=toc;
    time(2,2,i) =  t;
    npoints(2,2,i) = out_param.npoints;
    tic;
    [~, out_param] = funappx_g(f3,a,b,abstol,100,1000);
    t=toc;
    time(3,1,i) =   t;
    npoints(3,1,i) = out_param.npoints;
    tic;
    [~, out_param] = funappxglobal_g(f3,a,b,abstol,100,1000);
    t=toc;
    time(3,2,i) =   t;
    npoints(3,2,i) = out_param.npoints;
    tic;
    [~, out_param] = funappx_g(f4,a,b,abstol,100,1000);
    t=toc;
    time(4,1,i) = t;
    npoints(4,1,i) = out_param.npoints;
    tic;
    [~, out_param] = funappxglobal_g(f4,a,b,abstol,100,1000);
    t=toc;
    time(4,2,i) =   t;
    npoints(4,2,i) = out_param.npoints;
    tic;
end;
warning('on','MATLAB:funappxglobal_g:exceedbudget')
warning('on','MATLAB:funappx_g:peaky')
warning('on','MATLAB:funappxglobal_g:peaky')

timeratio = zeros(nrep,n);
npointsratio = zeros(nrep,n);
for i=1:nrep;
    for j=1:n;
        timeratio(i,j) = time(j,1, i)/time(j,2, i);
    end
end

for i=1:nrep;
    for j=1:n;
        npointsratio(i,j) = npoints(j,1, i)/npoints(j,2, i); 
    end
end

subplot(1,2,1); 
plot(1:nrep*n,sort(timeratio(:)),'blue',1:nrep*n,ones(nrep*n,1),'red')
title('time ratio');
subplot(1,2,2); 
plot(1:nrep*n,sort(npointsratio(:)),'blue',1:nrep*n,ones(nrep*n,1),'red')
title('npoints ratio');

display(' ')
display('   Test        # Points          Time Used')
display(' Function   Local    Global     Local    Global')
for i=1:n
    display(sprintf(['%8.0f %8.0f  %7.0f %10.6f  %10.6f'],...
        [i mean(npoints(i,1,:)) mean(npoints(i,2,:)) mean(time(i,1,:)) mean(time(i,2,:))])) 
end
