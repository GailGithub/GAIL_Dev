function funappxseashell (a, b, c, n, azimuth, elevation, res)
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
%   T. Davis & K. Sigmon, MATLAB Primer, 7th edition, CRC Press, 2005, pp. 80.
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

% Copyright 2006, Tim Davis, University of Florida

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
t = linspace(0, 2*pi, res) ;
[u,v] = meshgrid(t) ;

% define the surface
x = (a*(1-v/(2*pi)).*(1+cos(u)) + c) .* cos(n*v) ;
y = (a*(1-v/(2*pi)).*(1+cos(u)) + c) .* sin(n*v) ;
z = b*v/(2*pi) + a*(1-v/(2*pi)) .* sin(u) ;

% plot the surface
figure(1)
surf(x,y,z,y)
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

xapprox = zeros(res,res);
yapprox = zeros(res,res);
zapprox = zeros(res,res);
f1 = @(x,y) (a*(1-y/(2*pi)).*(1+cos(x))+c).*cos(n*y);
f2 = @(x,y) (a*(1-y/(2*pi)).*(1+cos(x))+c).*sin(n*y);
f3 = @(x,y) b*y/(2*pi)+a*(1-y/(2*pi)).*sin(x);
for i = 1:res;
    [xfappx,~]=funappx_g(@(x) f1(x,t(i)),0,2*pi);
    xapprox(i,:) = xfappx(t);
    [yfappx,~]=funappx_g(@(x) f2(x,t(i)),0,2*pi);
    yapprox(i,:) = yfappx(t);
    [zfappx,~]=funappx_g(@(x) f3(x,t(i)),0,2*pi);
    zapprox(i,:) = zfappx(t); 
end;

% plot the surface
figure(2)
surf(xapprox,yapprox,zapprox,yapprox)
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

errmat =sqrt((x-xapprox).^2+(y-yapprox).^2+(z-zapprox).^2);
figure(3);
surf(x,y,z,errmat);
colorbar;



