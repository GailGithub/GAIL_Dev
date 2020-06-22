function [fappx,npoints,errest] = funappx_g_gui(f,a,b,tol,nlo,nhi,varargin)
%funappx_g_gui Demonstrate numerical approximation of an univaraite function.
%   fappx = funappx_g_gui(f,a,b) shows the steps in approximating the function
%   f(x) from a to b by locally adaptive guaranteed method.
%   The color switches to green when the desired accuracy is obtained.
%
%   fappx = funappx_g_gui(f,a,b,tol) uses the given tolerance instead of 1.e-3
%   and returns an approximated function fappx%
%
%   [fappx, npoints] = funappx_g_gui(f,...) also gives the number points
%   needed of approximation.
%
%   Examples:
%  [fappx,npoints,errest] = funappx_g_gui(@(x) x.^2,-1,1,1e-2,10,20)
%  [fappx,npoints,errest] = funappx_g_gui(@(x) exp(-1000*(x-0.2).^2),0,1,1e-3,5,9)
%  [fappx,npoints,errest] = funappx_g_gui(@(x) exp(-1000*(x-0.2).^2),0,1,1e-6,10,20)
%   Flat function:
%  [fappx,npoints,errest] = funappx_g_gui(@(x) exp(-1./(x - 0.5).^2),0,1,1e-4,2,2)
%  [fappx,npoints,errest] = funappx_g_gui(@(x) sin(2*pi*x),0,1,1e-3,10,20)
%  Two local min:
%  [fappx,npoints] = funappx_g_gui(@(x) -5 * exp(-(10*(x - .3)).^2) - exp(-(10*(x - 0.75)).^2),0,1,1e-3,10,20)
%  [fappx,npoints] = ... 
%  Demo with funappxPenalty_g:
%  [fappx,npoints,errest] = funappx_g_gui(@(x) x.^2,-1,1,1e-2,10,20,'funappxPenalty_g')
shg
clf reset
MATLABVERSION = gail.matlab_version;
if isempty(varargin)
  algoname = 'funappx_g';
  algo = @(f,in_param) funappx_g(f,in_param);
else 
  algoname= varargin{1};
  algo = str2func(['@(f,in_param)', varargin{1},'(f,in_param)']);  
end
warning('off', ['GAIL:',algoname, ':exceediter']);


% Default tolerance
if nargin < 5
    nlo = 10;
    nhi = 1000;
end

if nargin < 4
    tol = 1.e-7;
    nlo = 10;
    nhi = 1000;
end

% Default function and interval.
if nargin < 3
    f = @(x) x.^2;
    a = 0;
    b = 1;
    tol = 1.e-3;
    nlo = 10;
    nhi = 1000;
end

% Initialization
fa = f(a);
fb = f(b);
k = 0;

% Scale the plot
h = b - a;
%ninit = 2*ceil(nhi*(nlo/nhi)^(1/(1+h)))+1;
ninit = nlo;
x = a:h/ninit:b;
x0 = x;
y = f(x);
maxy = max(y);
miny = min(y);
set(gcf,'doublebuffer','on','userdata',0)
MATLABBlue = [0, 0.447, 0.741];
MATLABGreen = [0.466,  0.674, 0.188];

hold on
p(1) = fill(a,fa,'k');
p(2) = fill(b,fb,'k');
hold off
s = (maxy - miny)/5;
axis([a b miny-s maxy+s])
% q(1) = uicontrol('string','step', ...
%     'units','normal','pos',[.65 .02 .08 .04], ...
%     'callback','set(gcf,''userdata'',1)');
% q(2) = uicontrol('string','auto', ...
%     'units','normal','pos',[.75 .02 .08 .04], ...
%     'callback','set(gcf,''userdata'',2)');
% q(3) = uicontrol('string','quit', ...
%     'units','normal','pos',[.85 .02 .08 .04], ...
%     'callback','set(gcf,''userdata'',3)');
q(1) = uicontrol('string','step', ...
    'units','normal','pos',[.75 .02 .08 .04], ...
    'callback','set(gcf,''userdata'',1)');
q(2) = uicontrol('string','auto', ...
     'units','normal','pos',[.85 .02 .08 .04], ...
    'callback','set(gcf,''userdata'',2)');
err = tol+1;

% Plot the input function
delta = 0.00001;
xx=a:delta:b; 


in_param.a = a; 
in_param.b = b; 
in_param.abstol = tol; 
in_param.nlo = nlo; 
in_param.nhi = nhi; 
in_param.output_x = 1;
in_param.ninit = ninit;
% tmpstr = strsplit(algoname,'_g');
% level = funmin_g(f,a,b,tol,nlo,nhi)-0.2;

while(max(err) > tol)
    if max(err) > tol;
        in_param.maxiter = k+1; 
        [~,out_param] = algo(f,in_param);
        err = out_param.errest;
        npoints = out_param.npoints;
        x = out_param.x;
        y = f(x);
        k = k + 1;
        p = flipud(get(gca,'children'));
        set(p(1),'xdata',x,'ydata',y)
        plot(x,y,'.','MarkerSize',20,'color',MATLABBlue);  hold on;
        plot(xx,f(xx),'color',MATLABBlue); 
        axis tight
        set(gca,'xtick', x0, 'xticklabel',x0);
        ax = gca;
        ax.XAxis.MinorTick = 'on';
        ax.XAxis.MinorTickValues = x;
        %plot(x,zeros(size(x)),'.','color',MATLABGreen); hold on;
        set(gca,'FontSize',16)
        title(['iter ', num2str(k), '; error bound is ' num2str(err),...
            '; number of points is ' num2str(npoints) ])
        %hTitle=title([tmpstr{1}, '\_g: error \(\approx\) ' sprintf('%0.2g',max(err)) ' in iter ' num2str(k)]);
        %set(hTitle,'FontSize',25,'Interpreter', 'latex')
        pause(.25)
        while get(gcf,'userdata') == 0
            pause(.25)
        end
        if get(gcf,'userdata') == 1
            set(gcf,'userdata',0)
        end
    else
        k = k + 1;
        break;
    end;
end

p = flipud(get(gca,'child'));
set(p(1),'xdata',x,'ydata',y)
set(gca,'xtick', x0, 'xticklabel',x0);
set(gca,'FontSize',16)
title(['iter ', num2str(k), '; error bound is ' num2str(err),...
    '; number of points is ' num2str(npoints) ])
%hTitle=title([tmpstr{1}, '\_g: error \(\approx\) '  sprintf('%0.2g',max(err)) ' in iter ' num2str(k)]);
%set(hTitle,'FontSize',25,'Interpreter', 'latex')
pause(.25)
while get(gcf,'userdata') == 0
    pause(.25)
end
if get(gcf,'userdata') == 1
    set(gcf,'userdata',0)
end
%npoints = index(end);
if MATLABVERSION >= 8.3
    fappx = griddedInterpolant(x,y,'linear');
else
    pp = interp1(x,y,'linear','pp');
    fappx =@(x) ppval(pp,x);
end;
errest = max(err);
%delete(p)
delete(q);
warning('on', ['GAIL:', algoname ,':exceediter']);


gail.save_eps('WorkoutFunappxOutput', [algoname, '_gui']);
%keyboard

% ---------------------------------------------------------
