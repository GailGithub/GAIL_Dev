function seashell (a, b, c, n, azimuth, elevation, res)
% SEASHELL draws a pretty Florida seashell, using a 3D parametric surface.
%
% Usage:
%
%   seashell (a, b, c, n, azimuth, elevation, res)
%   seashell ('spin') ;
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
%   If the azimuth is Inf, then the seashell view is spun dynamically.
%   Also try seashell ('spin') ;
%
% References:
%   T. Davis & K. Sigmon, MATLAB Primer, 8th edition, CRC Press, 2010.
%   von Seggern, CRC Standard Curves and Surfaces, 2nd edition, CRC Press,
%       1993, pp. 306-307.
%
% Example:
%   seashell ;          % draws the front cover of the MATLAB Primer
%   seashell (-0.5) ;   % draws the back cover
%   seashell (a,b,c,n,az,el,res) ;  % all options, defaults:
%                       % a=-0.2, b=0.5, c=0.1, n=2, az=-150, el=10, res=128
%
%   for a = -1:.1:1
%	seashell (a) ;
%	drawnow ;
%   end
%   for b = -1:.1:1
%	seashell (-.2, b) ;
%	drawnow
%   end
%
% See also SHELLGUI, SURF, VIEW, LINSPACE, MESHGRID, SHADING, LIGHTING,
%   LIGHTANGLE, COLORMAP, AXIS, MATERIAL, SIN, COS, PI.

% Copyright 2010, Tim Davis, University of Florida

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

% construct the res-by-res mesh
a1 = 0 ; b1 = 2 * pi;
t = a1: b1 / (res-1) : b1;
[u,v] = meshgrid(t) ;

% define the surface
f1 = @(x) cos(x);
f2 = @(x) sin(x);
x = (a*(1-v/(b1)).*(1+f1(u)) + c) .* f1(n*v) ;
y = (a*(1-v/(b1)).*(1+f1(u)) + c) .* f2(n*v) ;
z = (b/b1)*v + a*(1-v*(1/b1)) .* f2(u) ;

% plot the surface
% 7th Edition was surf(x,y,z,y)
surf(x,y,z,y-2*x)
shading interp

axis off
axis equal
colormap(hsv(1024))
material shiny
lighting gouraud
lightangle(80, -40)
lightangle(-90, 60)

% fix the view, or spin the seashell
if (isfinite (azimuth))
    view([azimuth elevation])
else
    for az = -180:10:180
        view ([az elevation])
        drawnow
    end
end

%% use funappx
% plot the surface

[cosappx, out1] = funappx_g( @(x) cos(x), a1, b1);
[sinappx, out2] = funappx_g( @(x) sin(x), a1, b1);
[cosappx2, out1] = funappx_g( @(x) cos(2*x), a1, 2*b1);
[sinappx2, out2] = funappx_g( @(x) sin(2*x), a1, 2*b1);
%u = repmat(t,res,1); v = u';
 
cosu  = repmat(cosappx(t),res,1);
sinu  = repmat(sinappx(t),res,1);
cosnv = repmat(cosappx2(t),res,1)'; 
sinnv = repmat(sinappx2(t),res,1)';
 
x1 = (a*(1-v/b1).*(1+cosu) + c) .* cosnv ;
y1 = (a*(1-v/b1).*(1+cosu) + c) .* sinnv ;
z1 = (b/b1)*v + a*(1-v/b1) .* sinu;

figure(2)
surf(x1,y1,z1,y1-2*x1)
shading interp

axis off
axis equal
colormap(hsv(1024))
material shiny
lighting gouraud
lightangle(80, -40)
lightangle(-90, 60)

% fix the view, or spin the seashell
if (isfinite (azimuth))
    view([azimuth elevation])
else
    for az = -180:10:180
        view ([az elevation])
        drawnow
    end
end

errmat = sqrt((x-x1).^2+(y-y1).^2+(z-z1).^2);
figure(5)
surf(x,y,z,errmat);
colorbar;
gail.save_eps('WorkoutFunappxOutput', 'Seashellsurfyerror');
