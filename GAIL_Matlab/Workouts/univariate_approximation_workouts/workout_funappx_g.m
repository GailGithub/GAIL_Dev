%comparison between funappx_g and funappxlocal_g
function [timelgratio,npointslgratio]=workout_funappx_g(nrep,abstol,nlo,nhi)
% user can choose absolut error tolerance, initial number of points, number
% of iteration or can use the following parameters
% nrep = 100; abstol = 1e-7; nlo = 100; nhi = 1000;
c = rand(nrep,1)*4;
n = 4;
npoints = zeros(n,2,nrep);
time = zeros(n,2,nrep);
warning('off','MATLAB:funappx_g:peaky')
warning('off','MATLAB:funappx_g:exceedbudget')
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
    [~, out_param] = funappx_g(f1,a,b,abstol,nlo,nhi);
    t=toc;
    time(1,1,i) = t;
    npoints(1,1,i) = out_param.npoints;
    tic;
    [~, out_param] = funappxglobal_g(f1,a,b,abstol,nlo,nhi);
    t=toc;
    time(1,2,i) = t;
    npoints(1,2,i) = out_param.npoints;
    tic;
    [~, out_param] = funappx_g(f2,a,b,abstol,nlo,nhi);
    t=toc;
    time(2,1,i) =  t;
    npoints(2,1,i) = out_param.npoints;
    tic;
    [~, out_param] = funappxglobal_g(f2,a,b,abstol,nlo,nhi);
    t=toc;
    time(2,2,i) =  t;
    npoints(2,2,i) = out_param.npoints;
    tic;
    [~, out_param] = funappx_g(f3,a,b,abstol,nlo,nhi);
    t=toc;
    time(3,1,i) =   t;
    npoints(3,1,i) = out_param.npoints;
    tic;
    [~, out_param] = funappxglobal_g(f3,a,b,abstol,nlo,nhi);
    t=toc;
    time(3,2,i) =   t;
    npoints(3,2,i) = out_param.npoints;
    tic;
    [~, out_param] = funappx_g(f4,a,b,abstol,nlo,nhi);
    t=toc;
    time(4,1,i) = t;
    npoints(4,1,i) = out_param.npoints;
    tic;
    [~, out_param] = funappxglobal_g(f4,a,b,abstol,nlo,nhi);
    t=toc;
    time(4,2,i) =   t;
    npoints(4,2,i) = out_param.npoints;
    tic;
end;
warning('on','MATLAB:funappxglobal_g:exceedbudget')
warning('on','MATLAB:funappx_g:peaky')
warning('on','MATLAB:funappx_g:exceedbudget')
warning('on','MATLAB:funappxglobal_g:peaky')

timeratio = zeros(nrep,n);
npointsratio = zeros(nrep,n);
for i=1:nrep;
    for j=1:n;
        timeratio(i,j) = time(j,1,i)/time(j,2, i);
    end
end

for i=1:nrep;
    for j=1:n;
        npointsratio(i,j) = npoints(j,1,i)/npoints(j,2, i); 
    end
end

timeratio = sort(timeratio(:));
npointsratio = sort(npointsratio(:));


%% Output the table
% To just re-display the output, load the .mat file and run this section
% only
display(' ')
display('   Test        # Points          Time Used')
display(' Function   Local    Global     Local    Global')
npointslgratio = zeros(1,n);
timelgratio = zeros(1,n);
for i=1:n
    display(sprintf('%8.0f %8.0f  %7.0f %10.6f  %10.6f',...
        [i mean(npoints(i,1,:)) mean(npoints(i,2,:)) mean(time(i,1,:)) mean(time(i,2,:))])) 
    npointslgratio(i) = mean(npoints(i,1,:))/mean(npoints(i,2,:));
    timelgratio(i) = mean(time(i,1,:))/mean(time(i,2,:));
end

%% Save Output

[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
if usejava('jvm')
figure
subplot(2,1,1);
plot(1:nrep*n,timeratio,'blue',1:nrep*n,ones(nrep*n,1),'red');
title('Comparison between funappx\_g and funappxglobal\_g')
ylabel('Time ratio of local/global')
xlabel('# of tests')
subplot(2,1,2);
plot(1:nrep*n,npointsratio,'blue',1:nrep*n,ones(nrep*n,1),'red');
ylabel('Points ratio of local/global')
xlabel('# of tests')
filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
    'WorkoutFunappxOutput',PATHNAMESEPARATOR,'WorkoutFunAppxTest-',...
    datestr(now,'dd-mmm-yyyy-HH-MM-SS'),'.eps');
print('-deps',filename)

filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
    'WorkoutFunappxOutput',PATHNAMESEPARATOR','WorkoutFunAppxTest-',...
    datestr(now,'dd-mmm-yyyy-HH-MM-SS'),'.mat');
save(filename, ...
    'npoints','time','c','timeratio','npointsratio','npointslgratio',...
    'timelgratio')
end
