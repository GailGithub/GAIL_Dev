function [fappx,npoints] = funappx_g_gui(f,a,b,tol,nlo,nhi)
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
%  [fappx,npoints]=funappx_g_gui(@(x) x.^2,-1,1,1e-2,10,20)
%  [fappx,npoints]=funappx_g_gui(@(x) exp(-1000*(x-0.2).^2),0,1,1e-2,10,20)
%   Flat function:
%  [fappx,npoints]=funappx_g_gui(@(x) exp(-1./(x - 0.5).^2),0,1,1e-5,10,20)
%  [fappx,npoints]=funappx_g_gui(@(x) sin(x),0,pi,1e-3,10,20)
%  [fappx,npoints]=funappx_g_gui(@(x) cos(x),0,pi,1e-3,10,20)
%  [fappx,npoints]=funappx_g_gui(@(x) sin(2*pi*x),0,1,1e-3,10,20)
%  Two local min:
%  [fappx,npoints]=funappx_g_gui(@(x) -5 * exp(-(10*(x - .3)).^2) - exp(-(10*(x - 0.75)).^2),0,1,1e-3,10,20)
%  
shg
clf reset
MATLABVERSION= gail.matlab_version;

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
index = [1 ninit];
% initialize nstar
%nstar = ninit - 2;
%nstar = floor(ninit/2);
% initialize error
err = tol + 1;
while(max(err) > tol)
    % length of each subinterval
    len = x(index(2:end))-x(index(1:end-1));
    reshapey = reshape(y(1:end-1),ninit - 1, (index(end) - 1)/(ninit -1));
    diffy = diff([reshapey;y(index(2:end))]);
    %approximate the weaker norm of input function at different subinterval
    %gn = (ninit-1)./len.*max(abs(diffy-repmat((y(index(2:end))-y(index(1:end-1)))/(ninit-1),ninit-1,1)),[],1);
    gn = (ninit-1)./len.*max(abs(bsxfun(@minus,diffy,(y(index(2:end))-y(index(1:end-1)))/(ninit-1))),[],1);
    %approximate the stronger norm of input function at different subinterval
    fn = (ninit-1)^2./(len.^2).*max(abs(diff(diffy)),[],1);
    %update cone condition every iteration
    ntemp=max(ceil(nhi*(nlo/nhi).^(1./(1+len))),3);
    nstar = ntemp -2;
    
    %     gn(gn<eps/2)=0;
    %     fn(fn<eps/2)=0;
    
    %find nstar not large enough then double it
    smallconeind = find(nstar.*(2*gn+fn.*len/(ninit-1)) <(fn.*len));
    nstar(smallconeind) = 2*nstar(smallconeind);
    
    err = nstar.*len.*gn./(4*(ninit-1).*(ninit-1-nstar));
    %check if error satisfy the error tolerance
    %     counterr = sum(err > abstol);
    %     if(length(x) + counterr *(ninit -1) > out_param.nmax)
    %         out_param.exit(1) = 1;
    %         warning('MATLAB:funappx_g:exceedbudget',' funappx_g attempted to exceed the cost budget. The answer may be unreliable.')
    %         break;
    %     end;
    if max(err) > tol;
        %flag sub interval error not satisfy error tolerance 1 in whbad
        whbad = err > tol;
        %add index for bad sub interval
        badind = find(whbad == 1);
        %flag sub interval error satisfy error tolerance 1 in whgood
        whgood = (whbad ==0);
        %add index for good sub interval
        goodind = find(whgood == 1);
        %find # of new sub intervals need to be added at each sub
        %interval
        badcumsum = cumsum(whbad);
        %pickup # of new sub intervals at bad intervals
        cumbad = badcumsum(badind);
        %generate new index of sub intervals splitted from bad intervals
        newindex = [badind + [0 cumbad(1:end-1)]; badind + cumbad];
        newindex = newindex(:)';
        %find the length of each sub interval
        h = len/2/(ninit-1);
        %reshape x without end point to a matrix of ninit-1 by # of intervals
        reshapex =  reshape(x(1:end-1),ninit -1,(index(end) - 1)/(ninit -1));
        %generate new points newx need to be added
        %newx = reshapex(:,badind) + repmat(h(badind),ninit-1,1);
        newx = bsxfun(@plus,reshapex(:,badind),h(badind));
        %compute value newy of newx
        newy = f(newx);
        %initialize a zero matrix of 2*(ninit-1) by # of bad sub intervals
        %to store all the points after splitting bad sub intervals
        badmatrix = zeros(2*(ninit-1),sum(whbad));
        %insert x at bad sub intervals in badmatrix as the row 1,
        %3,..., end-1
        badmatrix(1:2:end-1,:) = reshapex(:,badind);
        %insert newx at bad sub intervals in badmatrix as the row 2,
        %4,..., end
        badmatrix(2:2:end,:) = newx;
        %reshape badmatrix to the size of ninit -1 by 2*# of bad sub
        %intervals
        badmatreshape = reshape(badmatrix, ninit - 1, 2*sum(whbad));
        %initialize a matrix of ninit - 1 by # of sub intervals after
        %splitting bad sub intervals for x
        newreshapex = zeros(ninit - 1, 2*sum(whbad)+sum(whgood));
        %insert all the points after splitting bad sub intervals to correct
        %column
        newreshapex(:,newindex) = badmatreshape;
        %insert all the points on good sub intervals to correct column
        newreshapex(:,goodind + badcumsum(goodind)) = reshapex(:,goodind);
        %obtain all the points in vector x
        x = [newreshapex(:)' x(end)];
        %insert y at bad sub intervals in badmatrix as the row 1,
        %3,..., end-1
        badmatrix(1:2:end-1,:) = reshapey(:,badind);
        %insert newy at bad sub intervals in badmatrix as the row 2,
        %4,..., end
        badmatrix(2:2:end,:) = newy;
        %reshape badmatrix to the size of ninit -1 by 2*# of bad sub
        %intervals
        badmatreshape = reshape(badmatrix, ninit - 1, 2*sum(whbad));
        %initialize a matrix of ninit - 1 by # of sub intervals after
        %splitting bad sub intervals for y
        newreshapey = zeros(ninit - 1, 2*sum(whbad)+sum(whgood));
        %insert all the values after splitting bad sub intervals to correct
        %column
        newreshapey(:,newindex) = badmatreshape;
        %insert all the original y on good sub intervals to correct column
        newreshapey(:,goodind + badcumsum(goodind)) = reshapey(:,goodind);
        %obtain all the values in vector y
        y = [newreshapey(:)' y(end)];
        
        %generate error for new sub intervals
        %initialize a vertor of # of sub intervals after splitting
        newerr = zeros(1,2*sum(whbad)+sum(whgood));
        %use the same error for splitted bad interval
        baderr = [err(badind); err(badind)];
        %insert error after splitting bad sub intervals to correct
        %position
        newerr(newindex)=baderr(:)';
        newerr(goodind + badcumsum(goodind)) = err(goodind);
        %obtain error for all sub intervals
        err = newerr;
        %upadte index w.p.t x after splitting
        %update index of the original endpoints
        index(2:end) = index(2:end) + badcumsum*(ninit-1);
        %obtain the index of new endpoins after splitting
        %if one interval not splitted, will get the same index as in
        %previous line
        indexbeg = index(1:end-1) + whbad*(ninit-1);
        %combine two index together and emlinate duplicate indices
        indexnew = [index(1:end-1); indexbeg];
        indexnew = indexnew(:)';
        index = unique([indexnew index(end)]);
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
        title(['error estimation is ' num2str(max(err))])
        pause(.25)
        while get(gcf,'userdata') == 0
            pause(.25)
        end
        if get(gcf,'userdata') == 1
            set(gcf,'userdata',0)
        end
    else
        break;
    end;
end;
p = flipud(get(gca,'child'));
set(p(1),'xdata',x,'ydata',y)
set(gca,'xtick',x,'xticklabel',[]);
title(['error estimation is ' num2str(max(err))])
pause(.25)
while get(gcf,'userdata') == 0
    pause(.25)
end
if get(gcf,'userdata') == 1
    set(gcf,'userdata',0)
end
npoints = index(end);
if MATLABVERSION >= 8.3
    fappx = griddedInterpolant(x,y,'linear');
else
    pp = interp1(x,y,'linear','pp');
    fappx =@(x) ppval(pp,x);
end;
%delete(p)
delete(q);

% ---------------------------------------------------------
