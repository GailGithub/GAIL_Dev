function [npoints,errest,t,npointsglobal,errestglobal,tglobal]=funappx_convtest(f,a,b,varargin)
% a = -4.5960, b = 4.5960, c = 3.5960, f = @(x)exp(-1000*(x-c).^2)
% a = -4.4221, b = 4.4221, c = 3.4221, f = @(x)sin(c*pi*x)
%
% Compare funappx_g with funappxglobal_g:
%   [npoints,errest,t,npointsglobal,errestglobal,tglobal]=funappx_convtest(f,a,b)
%
% Compare funappxPenalty_g with funappxglobal_g:
%   [npoints,errest,t,npointsglobal,errestglobal,tglobal]=funappx_convtest(f,a,b,'funappxPenalty_g');
LatexInterpreter
tol = zeros(1,15);
errest = tol;
npoints = tol;
nmax = 500000000;
t = tol;
errestglobal = tol;
npointsglobal = tol;
tglobal = tol;
j=1;
%f = @(x) x.^2; 

if nargin > 1,
  in_param.a = a;
  in_param.b = b;
end

if nargin == 0,
  f = @(x) exp(-100*(x-sqrt(2)/2).^2);
  nmax = 1e7;
end

if isempty(varargin)
  algoname = 'funappx_g';
  algo = @(f,in_param) funappx_g(f,in_param);
else 
  algoname= varargin{1};
  algo = str2func(['@(f,in_param)', varargin{1},'(f,in_param)']);  
end

tmpstr = strsplit(algoname,'_g');
   
warning('off',['GAIL:',algoname,':peaky'])
warning('off',['GAIL:',algoname,':exceedbudget'])
warning('off','GAIL:funappxglobal_g:peaky')
warning('off','GAIL:funappxglobal_g:exceedbudget')

for i=-15:-1,
  tol(j) = 10^(i);
  in_param.abstol = 10^(i);
  in_param.nmax = nmax;
  tic,
  [~, out_param] = algo(f, in_param);
  t(j) = toc;
  errest(j) = out_param.errest;
  npoints(j) = out_param.npoints;
  
  tic,
  [~, out_param] = funappxglobal_g(f, in_param);
  tglobal(j) = toc;
  errestglobal(j) = out_param.errest;
  npointsglobal(j) = out_param.npoints;
  
  j=j+1;
end

warning('on','GAIL:funappxglobal_g:exceedbudget')
warning('on','GAIL:funappxglobal_g:peaky')
warning('on',['GAIL:',algoname,':peaky'])
warning('on',['GAIL:',algoname,':exceedbudget'])

[~,~,MATLABVERSION] = GAILstart(false); 
if usejava('jvm') || MATLABVERSION <= 7.12
    figure(1)
     
    semilogy(npoints, errest,'-o')
    title('Error estimate vs. computational cost')
    ylabel('error estimation')
    xlabel('Number of points')
    hold on
%     subplot(2,1,2)
%     semilogy(t, errest,'r--x')
%     ylabel('error estimation')
%     xlabel('time cost')
%     hold off
    
    %figure(2)
    semilogy(npointsglobal, errestglobal,'r--x')
    %title('Computational Cost of funappxglobal\_g VS error tolerance')
    ylabel('error estimation')
    xlabel('Number of points')
    legend([tmpstr{1}, '\_g'],'funappxglobal\_g')
   
    
%     subplot(2,1,2)
%     semilogy(tglobal, errestglobal,'r--x')
%     ylabel('error estimation')
%     xlabel('time cost')
%     legend('funappx\_g','funappxglobal\_g')
    
    hold off
    gail.save_eps('WorkoutFunappxOutput', ['Workout', algoname, 'ConvTest']);
end
end
