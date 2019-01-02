% funappx_g_test: comparison between funappx_g, funappxPenalty_g, and Chebfun
function [timeratio,npointsratio,matfilename] ...
   =funappx_g_test(nrep,abstol,varargin)
% user can choose absolute error tolerance, initial number of subintervals,
% number of iteration or can use the following parameters
% nrep = 100; abstol = 1e-6; ninit = 250;
% Compare funappx_g with funappxglobal_g and chebfun:
% [timeration,pointsratio]=funappx_g_test(nrep,abstol,ninit);
%
% Compare funappxPenalty_g with funappxglobal_g and chebfun:
% [timeratio,npointsratio]=funappx_g_test(nrep,abstol,ninit,'funappxPenalty_g');
rng default % for reproducibility
if nargin < 2
   abstol = 1e-6;
   if nargin < 1
      nrep = 1000;
   end
end
cc = rand(nrep,1);
c = cc*2; % number of simulations for each test function
n = 3; % number of test functions
m = 3; % number of methods
npoints = zeros(n,m,nrep);
time = zeros(n,m,nrep);
trueerrormat = zeros(n,m,nrep);
exceedmat  = zeros(n,m,nrep);
if isempty(varargin)
  algoname = 'funappx_g';
  algo = @(f,a,b,abstol) funappx_g(f,a,b,abstol);
elseif size(varargin) == 1
  algoname = 'funappx_g';
  algo = @(f,a,b,abstol) funappx_g(f,a,b,abstol,varargin{1});
else 
  algoname = varargin{2};
  algo = str2func(['@(f,a,b,abstol)', algoname,'(f,a,b,abstol)']);
end

algo2 = @(f,a,b,abstol) funappxglobal_g(f,a,b,abstol);
algo3 = @(f,a,b) chebfun(f, [a,b],'chebfuneps', abstol,'splitting','on');
methods = {algo, algo2, algo3};

warning('off',['GAIL:',algoname,':peaky'])
warning('off',['GAIL:',algoname,':exceedbudget'])
warning('off',['GAIL:',algoname,':fSmallerThanAbstol'])
warning('off','GAIL:funappxglobal_g:peaky')
warning('off','GAIL:funappxglobal_g:exceedbudget')

if ~exist('chebfun','file') 
   warning('Chebfun is not installed.')
   return
end
 
g1 = @(x,c) x.^4.*sin(c./((x==0)+x)); 
g2 = @(x,c) g1(x,c) + 10.*x.^2;
delta = .2; B = 1./(2*delta.^2);
g3 = @(x,cc) B*(4*delta.^2 + (x-cc).^2 + (x-cc-delta).*abs(x-cc-delta) ...
  - (x-cc+delta).*abs(x-cc+delta)).*(abs(x-cc) <= 2*delta);
%g4 = @(x,c) (x-c).^2;
%g5 = @(x,c) c*sin(c*pi*x);
%g6 = @(x,c) 10*exp(-1000*(x-c).^2);

a = zeros(1,n);
b = zeros(1,n);
a(1:3) = [-1,-1,-1];
b(1:3) = [1,1,1];
for i = 1:nrep
  gail.print_iterations(i, strcat(['Starting case i of ', int2str(nrep), ', i']), true);
  f1 = @(x) g1(x,c(i));
  f2 = @(x) g2(x,c(i));
  f3 = @(x) g3(x,cc(i)*0.6);
  %f4 = @(x) g4(x,c(i));
  %f5 = @(x) g5(x,c(i));
  %f6 = @(x) g6(x,c(i));
  fcns = {f1, f2, f3};
  %          f4 = @(x) 1/4*c(i)*exp(-2*x).*(c(i)-2*exp(x).*(-1 +...
  %              c(i)*cos(x) - c(i)*sin(x))+exp(2*x).*(c(i) + 2*cos(x)...
  %              - 2* sin(x) - c(i)*sin(2*x)));
  for j = 1:length(fcns)
    f = fcns{j};
    if j > 3
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
          tic, fappx = chebfun(f,[a(j),b(j)],'chebfuneps', abstol,'splitting','on'); t=toc;
          if ~isempty(lastwarn)
            exceedmat(j,k,i) = 1;
          end
        catch
            disp('Chebfun() error.')
            break
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
%       if trueerrormat(j,k,i)  > abstol
%         cf_chebfun(f,a(j),b(j), abstol);
%         keyboard
%       end
      if k==1
        exceedmat(j,k,i) = sum(out_param.exit)>0;
      elseif k==2
        exceedmat(j,k,i) = out_param.exceedbudget;
      end
    end
  end
end
warning('on','GAIL:funappxglobal_g:exceedbudget')
warning('on','GAIL:funappxglobal_g:peaky')
warning('on',['GAIL:',algoname,':peaky'])
warning('on',['GAIL:',algoname,':exceedbudget'])
warning('on',['GAIL:',algoname,':fSmallerThanAbstol'])

permuted_index = [3, 1:2];%, 4:length(fcns)];
% 
% d = 1.9;
% cc = d/4;
% f1 = @(x) g1(x,d);
% f2 = @(x) g2(x,d);
% f3 = @(x) g3(x,cc*0.6);
% %f4 = @(x) g4(x,d);
% f5 = @(x) g5(x,d);
% f6 = @(x) g6(x,d);
% fcns = {f1, f2, f3, f4, f5, f6};
% x = cell(length(fcns));
% y = cell(length(fcns));
% for i=1:length(fcns)
%     if i > 3,
%         b(i) = d+1;
%     end
%     x{i} = a(i):0.05:d+1;
%     y{i} = fcns{i}(x{i});
% end

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

%% Save Output
[~,~,MATLABVERSION] = GAILstart(false);
gail.InitializeDisplay
MATLABBlue = [0, 0.447, 0.741];
MATLABOrange = [0.85,  0.325, 0.098];
MATLABPurple = [0.494,  0.184, 0.556];
MATLABGreen = [0.466,  0.674, 0.188];
MATLABDkOrange = [0.85,  0.325, 0.098]*0.6;
MATLABLtOrange = 0.5*[0.85,  0.325, 0.098] + 0.5*[1 1 1];
%markers = {'--go', ':r*', '-.b.', '-g+', '--ro', '-.b'};
markers = {MATLABBlue, MATLABOrange, MATLABPurple, MATLABGreen, MATLABDkOrange,MATLABLtOrange};
if usejava('jvm') || MATLABVERSION <= 7.12
%   figure
%   hold on
%   for i = permuted_index
%     plot(x{i},y{i},'color',markers{i}); 
%   end
%   hold off
%  
%  legend('\(f_3\)','\(g_1\)','\(g_2\)','\(g_3\)','\(g_4\)','\(g_5\)','Location','NorthWest')
%  xlabel('x')
%  gail.save_eps('TraubPaperOutput', [algoname,'_testfun']);
  
  figure
  t = ((1:nrep*n) -1/2)/(nrep*n);
  for k = 2:m-1
    %subplot(1,m-1,k)
    semilogx(sorted_timeratio(k,:),t,'color',markers{k*2-1}); hold on
    semilogx(sorted_npointsratio(k,:),t,'color',markers{2*k});  hold on
    xlabel('Ratios'); 
    ylabel('Probability')
    %title([algoname, ' vs. ', func2str( methods{k+1})])
  end
  hold off
%   h=legend('{\tt funappx\_g} time / {\tt funappxglobal\_g} time',...
%            '{\tt funappx\_g} vs. {\tt funappxglobal\_g} points ratio',...
%            '{\tt funappx\_g} vs. Chebfun time ratio',...
%            '{\tt funappx\_g} vs. Chebfun points ratio',...
%            'Location','NorthWest');
%   set(h, 'Interpreter', 'latex')
%   legend BOXOFF 
  gail.save_eps('TraubPaperOutput', [algoname,'_test']);
end
matfilename = gail.save_mat('TraubPaperOutput', [algoname,'_test'], true, npoints, ...
  time, c, timeratio, npointsratio, nrep, n, m, ...
  sorted_timeratio, sorted_npointsratio, ...
  trueerrormat, exceedmat, permuted_index, abstol, ...
  algoname);
end


%% Sample printout
% # of replications = 1000
%    Test         Number of Points                    Time Used                          Success (%)                                  Failure (%)
%   Function   ----------------------------    -------------------------------     --------------------------------------   ----------------------------------------
%              Local      Global    Chebfun    Local       Global      Chebfun     Local        Global         Chebfun       Local         Global        Chebfun
%                                                                                  No Warn Warn No Warn Warn   No Warn Warn  No Warn Warn  No Warn Warn  No Warn Warn
%         3      5659     46439        116      0.0084       0.0137      0.0403     100      0    100      0      0        0      0      0      0      0    100      0
%         1      3690     26265         43      0.0091       0.0119      0.0104      98      0    100      0      3        0      2      0      0      0     97      0
%         2     11835     97106         22      0.0114       0.0246      0.0061     100      0    100      0      3        0      0      0      0      0     97      0
% IdleTimeout has been reached.