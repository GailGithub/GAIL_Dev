%funmin_g_test: comparison between funmin_g, fminbnd, and chebfun
function [timeratio,npointsratio,matfilename]=funmin_g_test(nrep,abstol,varargin)
% user can choose absolut error tolerance, initial number of points, number
% of iteration or can use the following parameters
% nrep = 100; abstol = 1e-6;
%
% Compare funmin_g, fminbnd, and chebfun:
% [timeratio,npointsratio]=funmin_g_test(nrep,abstol,'funmin_g');
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
  algoname = 'funmin_g';
  algo = @(f,a,b,abstol) funmin_g(f,a,b,abstol);
else
  algoname = varargin{1};
  algo = str2func(['@(f,a,b,abstol)', algoname,'(f,a,b,abstol)']);
end

algo2 = @(f,a,b,abstol) fminbnd(f,a,b,abstol);
algo3 = @(f,a,b) chebfun(f, [a,b],'chebfuneps', abstol,'splitting','on');
methods = {algo, algo2, algo3};

% warning('off',['GAIL:',algoname,':peaky'])
% warning('off',['GAIL:',algoname,':exceedbudget'])
% warning('off',['GAIL:',algoname,':fSmallerThanAbstol'])
% warning('off','GAIL:funminglobal_g:peaky')
% warning('off','GAIL:funminglobal_g:exceedbudget')

if ~exist('chebfun','file') 
   warning('Chebfun is not installed.')
   return
end

g1 = @(x,c) x.^4.*sin(c./((x==0)+x)); 
g2 = @(x,c) g1(x,c) + 10.*x.^2;
delta = .2; B = 1./(2*delta.^2);
g3 = @(x,cc) -B*(4*delta.^2 + (x-cc).^2 + (x-cc-delta).*abs(x-cc-delta) ...
  - (x-cc+delta).*abs(x-cc+delta)).*(abs(x-cc) <= 2*delta);
% g4 = @(x,c) (x-c).^2;
% g5 = @(x,c) c*sin(c*pi*x);
% g6 = @(x,c) -10*exp(-1000*(x-c).^2);

a = zeros(1,n);
b = zeros(1,n);
a(1:3) = [-1,-1,-1];
b(1:3) = [1,1,1];
for i = 1:nrep
  gail.print_iterations(i, strcat(['Starting case i of ', int2str(nrep), ', i']), true);
  f1 = @(x) g1(x,c(i));
  f2 = @(x) g2(x,c(i));
  f3 = @(x) g3(x,cc(i)*0.6);
%   f4 = @(x) g4(x,c(i));
%   f5 = @(x) g5(x,c(i));
%   f6 = @(x) g6(x,c(i));
  fcns = {f1, f2, f3};
  
  for j = 1:length(fcns)
    
    f = fcns{j};
    
    if j > 3
      b(j) = c(i)+1;
    end
    
    if j == 1
      exactxmin = -1;
      exactfmin = f(exactxmin);
    elseif j == 2
      exactfmin = 0;
    elseif j==3
      exactfmin = -1;
%     elseif j==4
%       exactxmin = c(i);
%       exactfmin = 0;
%     elseif j==5
%       kk = 0:floor(b(j)*c(i)-.5);
%       if isempty(kk)
%         exactfmin = 0;
%       else
%         stat_pts = (kk + 0.5)/c(i);
%         exactfmin = min(f([0,stat_pts,b(j)]));
%       end
%     elseif j==6
%       exactfmin = -10;
    end
    
    for k = 1:length(methods)
      if k == 1
        tic; [fmin, out_param] = methods{k}(f,a(j),b(j),abstol); t=toc;
        npoints(j,k,i) = out_param.npoints;
      elseif k == 2
        lastwarn('')
        tic, [~, fmin, exitflag, out] = fminbnd(@(x) f(x), a(j), b(j), ...
          optimset ('TolX', abstol, 'MaxFunEvals', out_param.nmax, 'MaxIter', out_param.maxiter) ); t=toc;
        if ~isempty(lastwarn) || exitflag ~= 1 || out.iterations > out_param.maxiter
          exceedmat(j,k,i) = 1;
        end
        npoints(j,k,i) = out.funcCount;
      elseif k == 3
        try
          lastwarn('')
          tic, chebf = chebfun(f,[a(j),b(j)],'chebfuneps', abstol,'splitting','on');
               fmin = min(chebf); 
          t=toc;
          if ~isempty(lastwarn)
            exceedmat(j,k,i) = 1;
          end
        catch
           disp('oops')
        end
        npoints(j,k,i) = length(chebf);
      end
      time(j,k,i) = t;
      
      trueerrormat(j,k,i) = max(abs(fmin-exactfmin));
%       if trueerrormat(j,k,i)  > abstol
%         cf_chebfun_min(f,a(j),b(j),abstol,exactfmin);
%         j, k
%         keyboard
%       end
      if k==1
        exceedmat(j,k,i) = sum(out_param.exitflag)>0;
      end
    end
  end
end
% warning('on','GAIL:funminglobal_g:exceedbudget')
% warning('on','GAIL:funminglobal_g:peaky')
% warning('on',['GAIL:',algoname,':peaky'])
% warning('on',['GAIL:',algoname,':exceedbudget'])
% warning('on',['GAIL:',algoname,':fSmallerThanAbstol'])

permuted_index = [3, 1:2];
% permuted_index = [3, 1:2, 4:length(fcns)];
% d = 1.9;
% cc = d/4;
% f1 = @(x) g1(x,d);
% f2 = @(x) g2(x,d);
% f3 = @(x) g3(x,cc*0.6);
% f4 = @(x) g4(x,d);
% f5 = @(x) g5(x,d);
% f6 = @(x) g6(x,d);
% fcns = {f1, f2, f3, f4, f5, f6};
% x = cell(length(fcns));
% y = cell(length(fcns));
% for i=1:length(fcns)
%   if i > 3,
%     b(i) = d+1;
%   end
%   x{i} = a(i):0.05:b(i);
%   y{i} = fcns{i}(x{i});
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
markers = {MATLABBlue, MATLABOrange, MATLABPurple, MATLABGreen, MATLABDkOrange,MATLABLtOrange};
if usejava('jvm') || MATLABVERSION <= 7.12
%   figure
%   hold on
%   for i = permuted_index
%     plot(x{i},y{i},'color',markers{i}); 
%   end
%   
%   legend('\(f_3\)','\(g_1\)','\(g_2\)','\(g_3\)','\(g_4\)','\(g_5\)','Location','NorthWest','Interpreter','latex')
%   xlabel('x')
%   gail.save_eps('TraubPaperOutput', [algoname,'_testfun']);
  
  figure
  t = ((1:nrep*n) -1/2)/(nrep*n);
  for k = 1:m-1
    semilogx(sorted_timeratio(k,:),t,'color',markers{k*2-1}); hold on
    semilogx(sorted_npointsratio(k,:),t,'color',markers{2*k});  hold on
    xlabel('Ratios'); 
    ylabel('Probability')
  end
%   hold off
%   h=legend('{\tt funmin\_g} vs. {\tt fminbnd} time ratio', '{\tt funmin\_g} vs. {\tt fminbnd} points ratio',...
%          '{\tt funmin\_g} vs. Chebfun time ratio', '{\tt funmin\_g} vs. Chebfun points ratio',...
%          'Location','NorthWest');
%   set(h, 'Interpreter', 'latex')   
%   legend BOXOFF 
  gail.save_eps('TraubPaperOutput', [algoname,'_test']);
end
matfilename = gail.save_mat('TraubPaperOutput', [algoname,'_test'], true, npoints, ...
  time, c, timeratio, npointsratio, nrep, n, m,...
  sorted_timeratio, sorted_npointsratio,...
  trueerrormat, exceedmat, permuted_index, abstol);
end

%% Sample printout
% # of replications = 1000
%    Test         Number of Points                    Time Used                          Success (%)                                  Failure (%)
%   Function   ----------------------------    -------------------------------     --------------------------------------   ----------------------------------------
%              funmin_g   fminbnd   Chebfun    funmin_g     fminbnd    Chebfun     funmin_g        fminbnd        Chebfun   funmin_g        fminbnd       Chebfun
%                                                                                  No Warn Warn No Warn Warn   No Warn Warn  No Warn Warn  No Warn Warn  No Warn Warn
%         3       274         8        116      0.0027       0.0019      0.0923     100      0    100      0     14        0      0      0      0      0     86      0 
%         1       230        22         43      0.0025       0.0025      0.0189     100      0     27      0     60        0      0      0     73      0     40      0 
%         2       273         9         22      0.0028       0.0021      0.0114     100      0    100      0     35        0      0      0      0      0     65      0 
