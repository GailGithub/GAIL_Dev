%traubpaper_funmin_g_test: comparison between funmin_g, fminbnd, and chebfun
function [timeratio,timelgratio,npointsratio,npointslgratio]=traubpaper_funmin_g_test(nrep,abstol,varargin)
% user can choose absolut error tolerance, initial number of points, number
% of iteration or can use the following parameters
% nrep = 10000; abstol = 1e-8;  
% 
% Compare funminNoPenalty_g, fminbnd, and chebfun:
% [timeratio,timelgratio,npointsratio,npointslgratio]=traubpaper_funmin_g_test(nrep,abstol,'funappxNoPenalty_g');

c = rand(nrep,1)*4; % number of simulations for each test function
n = 6; % number of test functions
m = 3; % number of methods
npoints = zeros(n,m,nrep);
time = zeros(n,m,nrep);
trueerrormat = zeros(n,m,nrep);
exceedmat  = zeros(n,m,nrep);
if isempty(varargin)
  algoname = 'funmin_g';
  algo = @(f,a,b,abstol) funmin_g(f,a,b,abstol);
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
                    lastwarn('') 
                    tic, fappx = chebfun(f,[a(j),b(j)]); t=toc;
                    if length(lastwarn) > 0
                      exceedmat(j,k,i) = 1;
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
 

timeratio = zeros(m-1,nrep,n);
npointsratio = zeros(m-1,nrep,n);
for i=1:nrep
    for j=1:n
        for k = 1:m-1
            timeratio(k,i,j) = time(j,1,i)/time(j,k+1,i);
        end
    end
end
for i=1:nrep;
    for j=1:n;
        for k = 1:m-1
            npointsratio(k,i,j) = npoints(j,1,i)/npoints(j,k+1,i);
        end
    end
end
sorted_timeratio = zeros(m-1,n*nrep);
sorted_npointsratio = zeros(m-1,n*nrep);
for k = 1:m-1
    sorted_timeratio(k,:) = sort(timeratio(k,:));
    sorted_npointsratio(k,:) = sort(npointsratio(k,:));
end

%% Output the table
% To just re-display the output, load the .mat file and run this section
% only
display(' ')
display('   Test         Number of Points                    Time Used                          Success (%)                                  Failure (%)')
display('  Function   ----------------------------    -------------------------------     --------------------------------------   ----------------------------------------');
display('             Local      Global    Chebfun    Local       Global      Chebfun     Local        Global         Chebfun       Local         Global        Chebfun')
display('                                                                                 No Warn Warn No Warn Warn   No Warn Warn  No Warn Warn  No Warn Warn  No Warn Warn')
npointslgratio = zeros(1,n);
timelgratio = zeros(1,n);

for i=1:n
    display(sprintf('%9.0f %9.0f %9.0f  %9.0f %11.3f  %11.3f %11.3f  %6.0f %6.0f %6.0f %6.0f %6.0f   %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f %6.0f',...
        [i mean(npoints(i,1,:)) mean(npoints(i,2,:)) mean(npoints(i,3,:))...
           mean(time(i,1,:)) mean(time(i,2,:)) mean(time(i,3,:))...
           100.0*sum(trueerrormat(i,1,:)<=abstol)/nrep 100.0*sum(trueerrormat(i,1,:)<=abstol & (exceedmat(i,1,:)))/nrep ...
           100.0*sum(trueerrormat(i,2,:)<=abstol)/nrep 100.0*sum(trueerrormat(i,2,:)<=abstol & (exceedmat(i,2,:)))/nrep ...
           100.0*sum(trueerrormat(i,3,:)<=abstol)/nrep 100.0*sum(trueerrormat(i,3,:)<=abstol & (exceedmat(i,3,:)))/nrep...
           100.0*sum(trueerrormat(i,1,:)>abstol)/nrep  100.0*sum(trueerrormat(i,1,:)>abstol & (exceedmat(i,1,:)))/nrep ...
           100.0*sum(trueerrormat(i,2,:)>abstol)/nrep  100.0*sum(trueerrormat(i,2,:)>abstol & (exceedmat(i,2,:)))/nrep ...
           100.0*sum(trueerrormat(i,3,:)>abstol)/nrep  100.0*sum(trueerrormat(i,3,:)>abstol & (exceedmat(i,3,:)))/nrep]))
    npointslgratio(i) = mean(npoints(i,1,:))/mean(npoints(i,2,:));
    timelgratio(i) = mean(time(i,1,:))/mean(time(i,2,:));
end

for k=1:m-1
    idx=find(timeratio(k,:,:)<1);
    max_idx_t = max(idx);
    timeratio(k,1:max_idx_t) = 1./timeratio(k,1:max_idx_t);
    
    idx=find(npointsratio(k,:,:)<1);
    max_idx_n = max(idx);
    npointsratio(k,1:max_idx_n) = 1.0 ./npointsratio(k,1:max_idx_n);
end

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
    gail.save_eps('TraubPaperOutput', ['traub_',algoname,'_testfun']);
    
    figure
    t =1:nrep*n;
    for k =1:m-1
        subplot(1,m-1,k)
        semilogy(t,sorted_timeratio(k,:),'r-',t,sorted_npointsratio(k,:),'b:');
        title([algoname, ' vs. ', func2str( methods{k+1})])
        legend('time ratio', 'points ratio');
    end
    
    gail.save_eps('TraubPaperOutput', ['traub_',algoname,'_test']);
end;
gail.save_mat('TraubPaperOutput', ['traub_',algoname,'_test'], true, npoints, ...
    time, c, timeratio, npointsratio, npointslgratio, timelgratio, ...
    sorted_timeratio, sorted_npointsratio);
end

%% Sample printout

%    Test         Number of Points                    Time Used                          Success (%)                                  Failure (%)
%   Function   ----------------------------    -------------------------------     --------------------------------------   ----------------------------------------
%              Local      Global    Chebfun    Local       Global      Chebfun     Local        Global         Chebfun       Local         Global        Chebfun
%                                                                                  No Warn Warn No Warn Warn   No Warn Warn  No Warn Warn  No Warn Warn  No Warn Warn
%         1     40021    331916       5364       0.100        0.058       1.461       0     98    100      0     26       32      0      2      0      0      0     42
%         2     42099    540010       5082       0.066        0.089       1.399      22     78    100      0     44       40      0      0      0      0      0     16
%         3      5948     97651          2       0.024        0.016       0.022      76     24    100      0     96        0      0      0      0      0      4      0
%         4     41926    608559          3       0.040        0.078       0.007       0    100    100      0    100        0      0      0      0      0      0      0
%         5    305060   1412518         39       0.482        0.196       0.013     100      0    100      0    100        0      0      0      0      0      0      0
%         6    112346   5315281        122       0.435        0.716       0.076       0    100    100      0     72        0      0      0      0      0     28      0

