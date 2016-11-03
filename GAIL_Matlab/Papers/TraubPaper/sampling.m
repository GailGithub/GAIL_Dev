InitializeWorkspaceDisplay
close all
set(0,'defaultaxesfontsize',48,'defaulttextfontsize',48, ... %make font larger
      'defaultLineLineWidth',7, ... %thick lines
      'defaultLineMarkerSize',48) %big dots
delta = 0.2;
B = 1./(2*delta.^2); 
c = -0.2;
f3 = @(x) -B*(4*delta.^2 + (x-c).^2 + (x-c-delta).*abs(x-c-delta) ...
     - (x-c+delta).*abs(x-c+delta)).*(abs(x-c) <= 2*delta); 
a = -1; b = 1; abstol = 1e-2;  
nlo = 10;
nhi = nlo;
hf=figure(1);
set(hf, 'Position', [1          96        2560        1113]);
%[fappx,out1]=funappx_g(f3,a,b,abstol,nlo,nhi)
% [fappx,npoints1,errest1] = funappx_g_gui(f3,a,b,abstol,nlo,nhi,...
%   'funappx_g')
[fappx,npoints1,errest1] = funappx_g_gui(f3,a,b,abstol,nlo,nhi)
xlabel('\(x\)')
ylabel('\(f_3(x)\)')
axis tight
gail.save_eps('TraubPaperOutput', 'sampling-funappxg');

hf2=figure(2);
set(hf2, 'Position', [1          96        2560        1113]);
%[fmin,out2]=funmin_g(f3,a,b,abstol,nlo,nhi)
[fmin,npoints2,errest2] = funmin_g_gui(f3,a,b,abstol,nlo,nhi,...
  'funmin_g')
xlabel('\(x\)')
ylabel('\(f_3(x)\)')
gail.save_eps('TraubPaperOutput', 'sampling-funming');