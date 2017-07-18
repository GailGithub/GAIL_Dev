% Working ones
% Example 1: 
f = @(x) prod(x,2); hyperbox = [zeros(1,2);ones(1,2)];
q = cubLattice_g(f,hyperbox,'uniform',1e-5,0,'transform','C1sin'); exactsol = 1/4;

w.f= @(x) prod(x,2);
w.absTol=1e-5
w.relTol=0 
w.domain = [zeros(1,2);ones(1,2)];
w.transform = 'C1sin';
cubLattice_gCLASS(w);

% Example 2: (Not working)
f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; hyperbox = [-inf(1,3);inf(1,3)];
q = cubLattice_g(f,hyperbox,'normal',1e-3,1e-3,'transform','C1sin','shift',2^(-25)*ones(1,3)); exactsol = 1;

w.f= @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2;
w.absTol=1e-3
w.relTol=1e-3
w.measure='normal';
w.domain = [-inf(1,3);inf(1,3)];
w.transform = 'C1sin';
w.shift=2^(-25)*ones(1,3);
cubLattice_gCLASS(w);

% Example 3:
f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-ones(1,2);2*ones(1,2)];
q = cubLattice_g(f,hyperbox,'uniform',1e-3,1e-2,'transform','C1'); exactsol = (sqrt(pi)/2*(erf(2)+erf(1)))^2;

w.f= @(x) exp(-x(:,1).^2-x(:,2).^2); 
w.absTol=1e-3
w.relTol=1e-2
w.domain = [-ones(1,2);2*ones(1,2)];
w.transform = 'C1';
cubLattice_gCLASS(w);

% Example 4: 
f = @(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); hyperbox = [-inf(1,1);inf(1,1)];
q = cubLattice_g(f,hyperbox,'normal',1e-4,1e-2,'transform','C1sin'); price = normcdf(0.05)*100 - 0.5*100*exp(-0.05^2/2);

w.f = @(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0);
w.domain = [-inf(1,1);inf(1,1)];
w.measure='normal'
w.absTol=1e-4
w.relTol=1e-2
w.transform=('C1sin')
cubLattice_gCLASS(w);

% Example 5: 
f = @(x) 8*prod(x,2); hyperbox = [zeros(1,5);ones(1,5)];
q = cubLattice_g(f,hyperbox,'uniform',1e-5,0); exactsol = 1/4;

w.f = @(x) 8*prod(x,2);
w.domain = [zeros(1,5);ones(1,5)];
w.measure='uniform'
w.absTol=1e-5
w.relTol=0
cubLattice_gCLASS(w);

% Example 6:
f = @(x) 3./(5-4*(cos(2*pi*x))); hyperbox = [0;1];
q = cubLattice_g(f,hyperbox,'uniform',1e-5,0,'transform','id'); exactsol = 1;

w.f= @(x) 3./(5-4*(cos(2*pi*x)));
w.absTol=1e-5
w.relTol=0
w.domain = [0;1];
w.transform = 'id';
cubLattice_gCLASS(w);

% Example 8: 
w.func = @(x)[10*x(:,1)-5*x(:,2).^2+1*x(:,3).^3, x(:,1), x(:,2).^2];
w.cv = [8,32/3]; hyperbox= [zeros(1,3);2*ones(1,3)];
[q, hi]= cubLattice_g(w,hyperbox,'uniform',1e-6,0); exactsol=128/3;

w.f= @(x)[10*x(:,1)-5*x(:,2).^2+1*x(:,3).^3, x(:,1), x(:,2).^2];
w.domain = [zeros(1,3);2*ones(1,3)];
w.trueMuCV=[8,32/3];
w.absTol=1e-6
w.relTol=0
w.measure = 'uniform';
[~, bye]=cubLattice_gCLASS(w);

% ------------------------------------
% Example 7:
f = @(x) [x(:,1).^2+x(:,2).^2]; hyperbox = [0,0,1];
[q, hi] = cubLattice_g(f,hyperbox,'uniform ball','abstol',1e-4,'reltol',0); exactsol = pi/2;

a.f = @(x) [x(:,1).^2+x(:,2).^2];
a.domain =[0,0,1];
a.measure = 'uniform ball';
a.abstol=1e-4;
a.reltol=0;
cubLattice_gCLASS(a);





