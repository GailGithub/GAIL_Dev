function funappxseashell(a, b, c, n, azimuth, elevation, res)
% funappxseashell draws a pretty Florida funappxseashell, using a 3D parametric surface.
%
% Usage:
%
%   funappxseashell (a, b, c, n, azimuth, elevation, res)
%   funappxseashell ('spin') ;
%
%   All arguments are optional.  The first four control the coefficients of
%   the parametric surface (u and v are the surface parameters):
%
%       x = (a*(1-v/(2*pi)).*(1+cos(u)) + c) .* cos(n*v) ;
%       y = (a*(1-v/(2*pi)).*(1+cos(u)) + c) .* sin(n*v) ;
%       z = b*v/(2*pi) + a*(1-v/(2*pi)) .* sin(u) ;
%
%   a,b: these determine how pointy or flat the shell is (informally...)
%   c: determines how much the shell overlaps with itself
%   n: the number of spirals in the shell
%
%   azimuth, elevation: determines the viewing angle (see the 'view' function)
%   res: the mesh size (res-by-res).  A larger number gives a smoother surface.
%
%   If the azimuth is Inf, then the funappxseashell view is spun dynamically.
%   Also try funappxseashell ('spin') ;
%
% References:
%   T. Davis & K. Sigmon, MATLAB Primer, 8th edition, CRC Press, 2010.
%   von Seggern, CRC Standard Curves and Surfaces, 2nd edition, CRC Press,
%       1993, pp. 306-307.
%
% Example:
%   funappxseashell ;          % draws the front cover of the MATLAB Primer
%   funappxseashell (-0.5) ;   % draws the back cover
%   funappxseashell (a,b,c,n,az,el,res) ;  % all options, defaults:
%                       % a=-0.2, b=0.5, c=0.1, n=2, az=-150, el=10, res=128
%
%   for a = -1:.1:1
%	  funappxseashell (a) ;
%	  drawnow ;
%   end
%   for b = -1:.1:1
%	  funappxseashell (-.2, b) ;
%	  drawnow
%   end
%
% See also SHELLGUI, SURF, VIEW, LINSPACE, MESHGRID, SHADING, LIGHTING,
%   LIGHTANGLE, COLORMAP, AXIS, MATERIAL, SIN, COS, PI.

% Copyright 2010, Tim Davis, University of Florida

set(0,'defaultaxesfontsize',13,'defaulttextfontsize',13) %make font larger
%set(0,'defaultLineLineWidth',3) %thick lines
%set(0,'defaultTextInterpreter','latex') %latex axis labels
%set(0,'defaultLineMarkerSize',40) %latex axis labels

% use default input parameters, if not present
if (nargin == 1 && ischar (a))
    in = -1 ;
else
    in = nargin ;
end
if (in < 1)
    a = -0.2 ;
end
if (in < 2)
    b = 0.5 ;
end
if (in < 3)
    c = 0.1 ;
end
if (in < 4)
    n = 2 ;
end
if (in < 5)
    azimuth = -150 ;
end
if (in < 6)
    elevation = 10 ;
end
if (in < 7)
    res = 128 ;
end
if (in == -1)
    azimuth = Inf ;
end

% sanity checks
if (a == 0)
    a = 0.01 ;
end
if (n <= 0)
    n = 0.1 ;
end

tic,
% construct the res-by-res mesh
t = linspace(0, 2*pi, res) ;
[u,v] = meshgrid(t) ;

% define the surface
x = (a*(1-v/(2*pi)).*(1+cos(u)) + c) .* cos(n*v) ;
y = (a*(1-v/(2*pi)).*(1+cos(u)) + c) .* sin(n*v) ;
z = b*v/(2*pi) + a*(1-v/(2*pi)) .* sin(u) ;
toc


%% use funappx
a1 = 0 ; b1 = 2 * pi;
[cosappx, out1] = funappxPenalty_g( @(x) cos(x), a1, 2*b1, 1e-1, 'nhi', res);
[sinappx, out2] = funappxPenalty_g( @(x) sin(x), a1, 2*b1, 1e-1, 'nhi', res);

tic, 
t0  = a1: b1 / (res-1) : 2*b1;
w = cosappx(t0);
w2 = sinappx(t0);
%v = repmat((1/(2*pi))*t(1:res), res,1)';
% x1 = (a*(1-v).*(repmat(1+w(1:res),res,1)) + c) .* repmat(w(1:2:end),res,1)' ;
% y1 = (a*(1-v).*(repmat(1+w(1:res),res,1)) + c) .* repmat(w2(1:2:end),res,1)' ;
% z1 = b*v + a*(1-v) .* repmat(w2(1:res),res,1) ;
vv = v/(2*pi);
x1 = bsxfun(@times,a*bsxfun(@times,1-vv,1+w(1:res))+c,w(1:2:end)');
y1 = bsxfun(@times,a*bsxfun(@times,1-vv,1+w(1:res))+c,w2(1:2:end)');
z1 = b*vv +a*bsxfun(@times,1-vv,w2(1:res));
toc


errest_cos = out1.errest
errest_sin = out2.errest

% plot the surface
% 7th Edition was surf(x,y,z,y)
figure;
surf(x,y,z,y)
shading interp

axis off
axis equal
colormap(hsv(1024))
material shiny
lighting gouraud
lightangle(80, -40)
lightangle(-90, 60)

% fix the view, or spin the funappxseashell
if (isfinite (azimuth))
    view([azimuth elevation])
else
    for az = -180:10:180
        view ([az elevation])
        drawnow
    end
end


gail.save_eps('WorkoutFunappxOutput', 'funappxseashell');

% plot the surface from funappx_G
% 7th Edition was surf(x,y,z,y)
figure(2)
surf(x1,y1,z1,y1)
shading interp

axis off
axis equal
colormap(hsv(1024))
material shiny
lighting gouraud
lightangle(80, -40)
lightangle(-90, 60)


% fix the view, or spin the funappxseashell
if (isfinite (azimuth))
    view([azimuth elevation])
else
    for az = -180:10:180
        view ([az elevation])
        drawnow
    end
end

gail.save_eps('WorkoutFunappxOutput', 'funappxseashell');

% plot the error
figure(3)
errmat = sqrt((x-x1).^2+(y-y1).^2+(z-z1).^2);
surf(x,y,z,errmat);
axis equal
colorbar;

if (isfinite (azimuth))
    view([azimuth elevation])
else
    for az = -180:10:180
        view ([az elevation])
        drawnow
    end
end

gail.save_eps('WorkoutFunappxOutput', 'Seashellsurfyerror');
