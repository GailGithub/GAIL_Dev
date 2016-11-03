function [fappxglobal,npoints,errorbound] = funappxglobal_g_gui(f,a,b,tol,nlo,nhi,nmax)
%funappx_g_gui  Demonstrate numerical approximation of an univaraite function.
%   fappx = FUNAPPX_G_GUI(f,a,b) shows the steps in approximating the function
%   f(x) from a to b by locally adaptive guaranteed method.
%   The color switches to green when the desired accuracy is obtained.
%
%   fappx = FUNAPPX_G_GUI(f,a,b,tol) uses the given tolerance instead of 1.e-3
%   and returns an approximated function fappx%
%
%   [fappx, npoints] = FUNAPPX_G_GUI(f,...) also gives the number points
%   needed of approximation.
%
%   Examples:
%  [fappxglobal,npoints,errorbound]=funappxglobal_g_gui(@(x) x.^2,-1,1,1e-2,10,20,10000000)
%  [fappxglobal,npoints,errorbound]=funappxglobal_g_gui(@(x) exp(-1000*(x-0.2).^2),0,1,1e-3,10,20,10000000)
%  [fappxglobal,npoints,errorbound]=funappxglobal_g_gui(@(x) exp(-1000*(x-0.2).^2),0,1,1e-6,10,20,10000000)
%   Flat function:
%  [fappxglobal,npoints,errorbound]=funappxglobal_g_gui(@(x) exp(-1./(x - 0.5).^2),0,1,1e-5,5,5,10000000)
%  [fappxglobal,npoints,errorbound]=funappxglobal_g_gui(@(x) sin(2*pi*x),0,1,1e-3,10,20,10000000)
%  Two local min:
%  [fappxglobal,npoints]=funappxglobal_g_gui(@(x) -5 * exp(-(10*(x - .3)).^2) - exp(-(10*(x - 0.75)).^2),0,1,1e-3,10,20)
%  
shg
clf reset
MATLABVERSION= gail.matlab_version;

% Default tolerance
% if nargin < 4
%     nmax = 10000000;
%     tol = 1.e-4;
%     nlo = 10;
%     nhi = 1000;
% end
% 
% % Default function and interval.
% if nargin < 3
%     f = @(x) x.^2;
%     a = 0;
%     b = 1;
%     tol = 1.e-3;
%     nlo = 10;
%     nhi = 1000;
% end


% Initialization
fa = f(a);
fb = f(b);

% Scale the plot
h = b - a;
ninit = max(ceil(nhi*(nlo/nhi)^(1/(1+h))),3);
x = a:h/(ninit-1):b;
y = f(x);
maxy = max(y);
miny = min(y);
set(gcf,'doublebuffer','on','userdata',0)
plot(x,y,'.','markersize',15);
hold on
p(1) = fill(a,fa,'k');
p(2) = fill(b,fb,'k');
hold off
s = (maxy - miny)/5;
axis([a b miny-s maxy+s])
q(1) = uicontrol('string','step', ...
    'units','normal','pos',[.65 .02 .08 .04], ...
    'callback','set(gcf,''userdata'',1)');
q(2) = uicontrol('string','auto', ...
    'units','normal','pos',[.75 .02 .08 .04], ...
    'callback','set(gcf,''userdata'',2)');
q(3) = uicontrol('string','quit', ...
    'units','normal','pos',[.85 .02 .08 .04], ...
    'callback','set(gcf,''userdata'',3)');
% index = [1 ninit];
% initialize nstar
nstar = ninit - 2;
n = ninit;
tauchange = 0;
flag = 0;
errorbound =1;
while n < nmax;
    
    if(flag==0)
        x = a:h/(n-1):b;
        y = f(x);        
    else
        xnew = repmat(x(1:end-1),m-1,1)+repmat((1:m-1)'*h/(n-1),1,(n-1)/m);
        ynew = f(xnew);
        xnew = [x(1:end-1); xnew];
        x = [xnew(:); x(end)]';
        ynew = [y(1:end-1); ynew];
        y = [ynew(:); y(end)]';
    end;
    %        tt = length(index)-1;
        p = flipud(get(gca,'children'));
        %         p = flipud(allchild(gca));
        set(p(1),'xdata',x,'ydata',y)
        %         for i =1:tt;
        %             if err(i) < tol;
        %                 color = [0 .5 0];
        %             else
        %                 color = [.6 .6 .6];
        %             end
        %             left = index(i);
        %             right = index(i+1);
        %             u = [x(left) x(left:right) x(right)];
        %             v = [0 y(left:right) 0];
        %             set(p(i+1),'xdata',u,'ydata',v,'facecolor',color)
        %         end;
        set(gca,'xtick',x,'xticklabel',[]);
        title(['error estimation is ' num2str(max(errorbound))])
        pause(.25)
        while get(gcf,'userdata') == 0
            pause(.25)
        end
        if get(gcf,'userdata') == 1
            set(gcf,'userdata',0)
        end
    diff_y = diff(y);
    %approximate the weaker norm of input function
    gn = (n-1)/h*max(abs(diff_y-(y(n)-y(1))/(n-1)));
    %approximate the stronger norm of input function
    fn = (n-1)^2/h^2*max(abs(diff(diff_y)));
    errorbound = fn*h^2/(8*(n-1)^2);
    % Stage 2: satisfy necessary condition
    if nstar*(2*gn+fn*h/(n-1)) >= fn*h;
        % Stage 3: check for convergence
        errbound = 4*tol*(n-1)*(n-1-nstar)...
            /nstar/h;
        % satisfy convergence
        if errbound >= gn;
             break;
        end;
        % otherwise increase number of points
        m = max(ceil(1/(n-1)*sqrt(gn*nstar*...
            h/4/tol)),2);
        n = m*(n-1)+1;
        flag = 1;
        % Stage2: do not satisfy necessary condition
    else
        % increase tau
        nstar = fn/(2*gn/h+fn/(n-1));
        % change tau change flag
        tauchange = 1;
        % check if number of points large enough
        if n >= nstar+2;
            % true, go to Stage 3
            errbound = 4*tol*(n-1)*(n-1-nstar)...
                /nstar/h;
            if errbound >= gn;
                break;
            end;
            m = max(ceil(1/(n-1)*sqrt(gn*nstar*...
                h/4/tol)),2);
            n = m*(n-1)+1;
            flag = 1;
        else
            % otherwise increase number of points, go to Stage 1
            n = 2 + ceil(nstar);
        end;
    end;
end;

if tauchange == 1;
    warning('GAIL:funappxglobal_g:peaky','This function is peaky relative to nlo and nhi. You may wish to increase nlo and nhi for similar functions.')
end;
errorbound = fn*h^2/(8*(n-1)^2);
p = flipud(get(gca,'child'));
set(p(1),'xdata',x,'ydata',y)
set(gca,'xtick',x,'xticklabel',[]);
title(['error estimation is ' num2str(max(errorbound))])
pause(.25)
while get(gcf,'userdata') == 0
    pause(.25)
end
if get(gcf,'userdata') == 1
    set(gcf,'userdata',0)
end
npoints = n;
if MATLABVERSION >= 8.3
    fappxglobal = griddedInterpolant(x,y,'linear');
else
    pp = interp1(x,y,'linear','pp');
    fappxglobal =@(x) ppval(pp,x);
end;
%delete(p)
delete(q);

% ---------------------------------------------------------
