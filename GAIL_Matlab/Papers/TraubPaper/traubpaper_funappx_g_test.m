% traubpaper_funappx_g_test: comparison between funappx_g and funappxlocal_g
function [timeratio,timelgratio,npointsratio,npointslgratio]=traubpaper_funappx_g_test(nrep,abstol,varargin)
% user can choose absolut error tolerance, initial number of points, number
% of iteration or can use the following parameters
% nrep = 100; abstol = 1e-6;
%
% Compare funappxNoPenalty_g with funappxglobal_g and chebfun:
% [timeratio,timelgratio,npointsratio,npointslgratio]=traubpaper_funappx_g_test(nrep,abstol,'funappxNoPenalty_g');
rng default % for reproducibility
set(0,'defaultaxesfontsize',14,'defaulttextfontsize',14, ... %make font larger
  'defaultLineLineWidth',2); %thick lines
%  'defaultLineMarkerSize',8)
cc = rand(nrep,1);
c = cc*2; % number of simulations for each test function
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
algo3 = @(f,a,b) chebfun(f, [a,b],'splitting','on');
methods = {algo, algo2, algo3};

warning('off',['GAIL:',algoname,':peaky'])
warning('off',['GAIL:',algoname,':exceedbudget'])
warning('off',['GAIL:',algoname,':fSmallerThanAbstol'])
warning('off','GAIL:funappxglobal_g:peaky')
warning('off','GAIL:funappxglobal_g:exceedbudget')

g1 = @(x,c) x.^4 .* sin(c./x);
g2 = @(x,c) g1(x,c) + c.*x.^2;
delta = .2; B = 1./(2*delta.^2);
g3 = @(x,cc) B*(4*delta.^2 + (x-cc).^2 + (x-cc-delta).*abs(x-cc-delta) ...
  - (x-cc+delta).*abs(x-cc+delta)).*(abs(x-cc) <= 2*delta);
g4 = @(x,c) (x-c).^2;
g5 = @(x,c) c*sin(c*pi*x);
g6 = @(x,c) 10*exp(-1000*(x-c).^2);

a = zeros(1,n);
b = zeros(1,n);
a(1:3) = [-1,-1,-1];
b(1:3) = [1,1,1];
for i = 1:nrep
  f1 = @(x) g1(x,c(i));
  f2 = @(x) g2(x,c(i));
  f3 = @(x) g3(x,cc(i)*0.6);
  f4 = @(x) g4(x,c(i));
  f5 = @(x) g5(x,c(i));
  f6 = @(x) g6(x,c(i));
  fcns = {f1, f2, f3, f4, f5, f6};
  %          f4 = @(x) 1/4*c(i)*exp(-2*x).*(c(i)-2*exp(x).*(-1 +...
  %              c(i)*cos(x) - c(i)*sin(x))+exp(2*x).*(c(i) + 2*cos(x)...
  %              - 2* sin(x) - c(i)*sin(2*x)));
  
  
  for j = 1:length(fcns)
    f = fcns{j};
    if j > 3,
      b(j) = c(i)+1;
    end
    xx = a(j):0.000001:b(j);%rand(1000000,1)*(b-a)+a;
    exactyy = f(xx);
    for k = 1:length(methods)
      clear fappx out_param yy t
      if k <= 2
        tic; [fappx, out_param] = methods{k}(f,a(j),b(j),abstol); t=toc;
        npoints(j,k,i) = out_param.npoints;
      elseif k == 3
        try
          lastwarn('')
          tic, fappx = chebfun(f,[a(j),b(j)],'splitting','on'); t=toc;
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

permuted_index = [3, 1:2, 4:length(fcns)];
d = 1.9;
cc = d/4;
f1 = @(x) g1(x,d);
f2 = @(x) g2(x,d);
f3 = @(x) g3(x,cc*0.6);
f4 = @(x) g4(x,d);
f5 = @(x) g5(x,d);
f6 = @(x) g6(x,d);
fcns = {f1, f2, f3, f4, f5, f6};
x = cell(length(fcns));
y = cell(length(fcns));
for i=1:length(fcns)
    if i > 3,
        b(i) = d+1;
    end
    x{i} = a(i):0.05:d+1;
    y{i} = fcns{i}(x{i});
end

timeratio = zeros(m-1,nrep,n);
npointsratio = zeros(m-1,nrep,n);
for i=1:nrep  % each test function
  for j=1:n % jth family of test functions
    for k = 1:m-1
      timeratio(k,i,j) = time(j,1,i)/time(j,k+1,i); % first method compared to k+1 th method
    end
  end
end
for i=1:nrep  % each test function
  for j=1:n % jth family of test functions
    for k = 1:m-1
      npointsratio(k,i,j) = npoints(j,1,i)/npoints(j,k+1,i); % first method compared to k+1 th method
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

for i = permuted_index
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
markers = {'--go', ':r*', '-.b.', '-g+', '--ro', '-.b'};
if usejava('jvm') || MATLABVERSION <= 7.12
  figure
  for i = permuted_index
    plot(x{i},y{i}, markers{i}); hold on
  end
  
  legend('f3','g1','g2','g3','g4','g5','Location','NorthWest')
  gail.save_eps('TraubPaperOutput', ['traub_',algoname,'_testfun']);
  
  figure
  t =1:nrep*n;
  for k =1:m-1
    subplot(1,m-1,k)
    semilogy(t,sorted_timeratio(k,:),'r-',t,sorted_npointsratio(k,:),'b:');
    %title([algoname, ' vs. ', func2str( methods{k+1})])
    legend('time ratio', 'points ratio','Location','NorthWest');
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
%         3      2904     46439         11       0.009        0.015       0.183     100      0    100      0     98        0      0      0      0      0      2      0
%         1      2924     27380       5249       0.015        0.011       2.092     100      0    100      0    100       71      0      0      0      0      0      0
%         2      3074     38644       1489       0.009        0.013       0.908     100      0    100      0    100        0      0      0      0      0      0      0
%         4      3493     26524          3       0.010        0.009       0.009     100      0    100      0    100        0      0      0      0      0      0      0
%         5      6621     37102         22       0.013        0.007       0.010     100      0    100      0    100        0      0      0      0      0      0      0
%         6     11398    345861        151       0.014        0.042       0.132     100      0    100      0    100        0      0      0      0      0      0      0
