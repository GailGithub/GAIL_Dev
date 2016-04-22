%comparison between funappx_g and funappxlocal_g
function [timeratio,timelgratio,npointsratio,npointslgratio]=workout_funappx_g(nrep,abstol,varargin)
% user can choose absolut error tolerance, initial number of points, number
% of iteration or can use the following parameters
% nrep = 100; abstol = 1e-13; nlo = 100; nhi = 1000;
% 
% Compare funappx_g with funappxglobal_g:
% [timeratio,timelgratio,npointsratio,npointslgratio]=workout_funappx_g(nrep,abstol);
%
% Compare funappxNoPenalty_g with funappxglobal_g:
% [timeratio,timelgratio,npointsratio,npointslgratio]=workout_funappx_g(nrep,abstol,'funappxNoPenalty_g');

c = rand(nrep,1)*4;
n = 6; % number of test functions
m = 3; % number of methods
npoints = zeros(n,m,nrep);
time = zeros(n,m,nrep);
trueerrormat = zeros(n,m,nrep);
exceedmat  = zeros(n,m,nrep);
if isempty(varargin)
  algoname = 'funappx_g';
  algo = @(f,a,b,abstol) funappx_g(f,a,b,abstol);
else 
  algoname = varargin{1};
  algo = str2func(['@(f,a,b,abstol)', algoname,'(f,a,b,abstol)']);  
end

algo2 = @(f,a,b,abstol) funappxglobal_g(f,a,b,abstol);
algo3 = @(f,a,b) chebfun(f, [a,b]);
methods = {algo, algo2, algo3};

warning('off',['GAIL:',algoname,':peaky'])
warning('off',['GAIL:',algoname,':exceedbudget'])
warning('off',['GAIL:',algoname,':fSmallerThanAbstol'])
warning('off','GAIL:funappxglobal_g:peaky')
warning('off','GAIL:funappxglobal_g:exceedbudget')

a = zeros(1,n);
b = zeros(1,n);
a(1:3) = [-1,-1,-1];
b(1:3) = [1,1,1];
for i = 1:nrep
    f1 = @(x) x.^4 .* sin(c(i)./x);
    f2 = @(x) f1(x) + c(i).*x.^2;
    delta = .2; B = 1./(2*delta.^2); cc = -c(i);
    f3 = @(x) B*(4*delta.^2 + (x-cc).^2 + (x-cc-delta).*abs(x-cc-delta) ...
        - (x-cc+delta).*abs(x-cc+delta)).*(abs(x-cc) <= 2*delta);
    f4 = @(x) (x-c(i)).^2;
    f5 = @(x) c(i)*sin(c(i)*pi*x);
    f6 = @(x) 10*exp(-1000*(x-c(i)).^2);
    %          f4 = @(x) 1/4*c(i)*exp(-2*x).*(c(i)-2*exp(x).*(-1 +...
    %              c(i)*cos(x) - c(i)*sin(x))+exp(2*x).*(c(i) + 2*cos(x)...
    %              - 2* sin(x) - c(i)*sin(2*x)));
    fcns = {f1, f2, f3, f4, f5, f6};

    for j = 1:length(fcns)
        f = fcns{j};
        if j > 3,     
            b(j) = c(i)+1;
        end
        xx = a(j):0.000001:b(j);%rand(1000000,1)*(b-a)+a;
        exactyy = f(xx);
        for k = 1:length(methods)
            if k <=2
                tic; [fappx, out_param] = methods{k}(f,a(j),b(j),abstol); t=toc;
                npoints(j,k,i) = out_param.npoints;
            elseif k==3 
                try
                    splitting on
                    tic, fappx = chebfun(f,[a(j),b(j)]); t=toc;
                    if length(lastwarn) > 0
                      exceedmat(j,k,i) = 1;
                      lastwarn('') 
                    end
                end
                npoints(j,k,i) = length(fappx);
            end
            time(j,k,i) = t;
            
            if k ==1 || k==3
                yy = fappx(xx);
            elseif k ==2
                yy = ppval(fappx,xx);  
            end
        
            trueerrormat(j,k,i) = max(abs(yy-exactyy));
%             if trueerrormat(j,k,i)  > abstol
%                cf_chebfun(f,a,b);
%                keyboard
%             end
            if k==1
                exceedmat(j,k,i) = sum(out_param.exitflag)>0;
            elseif k==2
                exceedmat(j,k,i) = out_param.exceedbudget;
            end
        end
    end
end;
warning('on','GAIL:funappxglobal_g:exceedbudget')
warning('on','GAIL:funappxglobal_g:peaky')
warning('on',['GAIL:',algoname,':peaky'])
warning('on',['GAIL:',algoname,':exceedbudget'])
warning('on',['GAIL:',algoname,':fSmallerThanAbstol'])

cc = 2.5;
x = cell(length(fcns));
y = cell(length(fcns));
for i=1:length(fcns)
   x{i}= a(i):0.05:b(i);
   y{i} = fcns{i}(x{i});
end
 

timeratio = zeros(nrep,n);
npointsratio = zeros(nrep,n);
for i=1:nrep
    for j=1:n
          timeratio(i,j) = time(j,1,i)/time(j,2,i);
    end
end

for i=1:nrep;
    for j=1:n;
        npointsratio(i,j) = npoints(j,1,i)/npoints(j,2,i); 
    end
end

timeratio = sort(timeratio(:));
npointsratio = sort(npointsratio(:));

%% Output the table
% To just re-display the output, load the .mat file and run this section
% only
display(' ')
display('   Test         Number of Points                    Time Used                          Success                                      Failure')
display('  Function   ----------------------------    -------------------------------     --------------------------------------   ----------------------------------------');
display('             Local      Global    Chebfun    Local       Global      Chebfun     Local        Global         Chebfun       Local         Global        Chebfun')
display('                                                                                 No Warn Warn No Warn Warn   No Warn Warn  No Warn Warn  No Warn Warn  No Warn Warn')
npointslgratio = zeros(1,n);
timelgratio = zeros(1,n);

for i=1:n
    display(sprintf('%9.0f %9.0f %9.0f  %9.0f %11.3f  %11.3f %11.3f  %6.0f %6.0f %6.0f %6.0f %6.0f   %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f',...
        [i mean(npoints(i,1,:)) mean(npoints(i,2,:)) mean(npoints(i,3,:))...
           mean(time(i,1,:)) mean(time(i,2,:)) mean(time(i,3,:))...
           100.0*sum(trueerrormat(i,1,:)<=abstol & (~exceedmat(i,1,:)))/nrep 100.0*sum(trueerrormat(i,1,:)<=abstol & (exceedmat(i,1,:)))/nrep ...
           100.0*sum(trueerrormat(i,2,:)<=abstol & (~exceedmat(i,2,:)))/nrep 100.0*sum(trueerrormat(i,2,:)<=abstol & (exceedmat(i,2,:)))/nrep ...
           100.0*sum(trueerrormat(i,3,:)<=abstol & (~exceedmat(i,3,:)))/nrep 100.0*sum(trueerrormat(i,3,:)<=abstol & (exceedmat(i,3,:)))/nrep...
           100.0*sum(trueerrormat(i,1,:)>abstol & (~exceedmat(i,1,:)))/nrep  100.0*sum(trueerrormat(i,1,:)>abstol & (exceedmat(i,1,:)))/nrep ...
           100.0*sum(trueerrormat(i,2,:)>abstol & (~exceedmat(i,2,:)))/nrep  100.0*sum(trueerrormat(i,2,:)>abstol & (exceedmat(i,2,:)))/nrep ...
           100.0*sum(trueerrormat(i,3,:)>abstol & (~exceedmat(i,3,:)))/nrep  100.0*sum(trueerrormat(i,3,:)>abstol & (exceedmat(i,3,:)))/nrep]))
    npointslgratio(i) = mean(npoints(i,1,:))/mean(npoints(i,2,:));
    timelgratio(i) = mean(time(i,1,:))/mean(time(i,2,:));
end

idx=find(timeratio<1);
max_idx_t = max(idx);
timeratio(1:max_idx_t) = 1./timeratio(1:max_idx_t);
idx=find(npointsratio<1);
max_idx_n = max(idx);
npointsratio(1:max_idx_n) = 1.0 ./npointsratio(1:max_idx_n);
 

%% Output the table
% To just re-display the output, load the .mat file and run this section
% only

%% Save Output
[~,~,MATLABVERSION] = GAILstart(false);
markers = {'--go', ':r*', '-.b*', '-g+', '--ro', '-.b'};
if usejava('jvm') || MATLABVERSION <= 7.12
    figure
    for i=1:length(fcns)
      plot(x{i},y{i}, markers{i}); hold on
    end

    legend('f1','f2','f3', 'f4', 'f5', 'f6')
    gail.save_eps('WorkoutFunappxOutput', 'testfun');
    
    figure
    t =1:nrep*n;
    plot(t,timeratio,'r',t,npointsratio,'b:');
    legend('time ratio','points ratio');
    gail.save_eps('WorkoutFunappxOutput', ['Workout',algoname,'Test']);
end;
gail.save_mat('WorkoutFunappxOutput', ['Workout',algoname,'Test'], true, npoints,time,...
    c,timeratio,npointsratio,npointslgratio,timelgratio);

end

%% If funappx_g is used:
% % Sample output for nrep=1000; abstol = 1e-7; nlo = 100; nhi = 1000;
% %    Test      Number of Points       Time Used
% %  Function   Local      Global     Local    Global
% %         1     72452     200725   0.0169790    0.0351434
% %         2    489265     401931   0.1527540    0.0790056
% %         3    225006    2091097   0.0763022    0.4285055

% % 
% % timelgratio =
% % 
% %        0.4831    1.9335    0.1781
% % 
% % 
% % npointslgratio =
% % 
% %        0.3610    1.2173    0.1076
%
%% If funappxNoPenalty_g is used:
% 
%    Test      Number of Points       Time Used
%  Function   Local      Global     Local    Global
%         1      7781     194166   0.0045961    0.0379585
%         2     53917     425327   0.0223668    0.0496360
%         3     22875    2092566   0.0158381    0.2506039
%
% timelgratio =
%
%     0.1211    0.4506    0.0632
%
% npointslgratio =
% 
%     0.0401    0.1268    0.0109

