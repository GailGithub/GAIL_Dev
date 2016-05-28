InitializeWorkspaceDisplay
close all
delta = 0.1;
B = 1./(2*delta.^2); 
c = -0.2;
f3 = @(x) -B*(4*delta.^2 + (x-c).^2 + (x-c-delta).*abs(x-c-delta) ...
     - (x-c+delta).*abs(x-c+delta)).*(abs(x-c) <= 2*delta); 
a = -1; b = 1; abstol = 11e-2;  
nlo = 10;
nhi = nlo;
figure(1)
%[fappx,out1]=funappxNoPenalty_g(f3,a,b,abstol,nlo,nhi)
[fappx,npoints1,errest1] = funappx_g_gui(f3,a,b,abstol,nlo,nhi,...
  'funappxNoPenalty_g')
%axis([a,b,-1.4,0.1])
xlabel('\(x\)')
ylabel('\(f_3\)')
gail.save_eps('TraubPaperOutput', 'sampling-funappxg');

figure(2)
%[fmin,out2]=funminNoPenalty_g(f3,a,b,abstol,nlo,nhi)
[fmin,npoints2,errest2] = funmin_g_gui(f3,a,b,abstol,nlo,nhi,...
  'funminNoPenalty_g')
%axis([a,b,-1.4,0.1])
gail.save_eps('TraubPaperOutput', 'sampling-funming');