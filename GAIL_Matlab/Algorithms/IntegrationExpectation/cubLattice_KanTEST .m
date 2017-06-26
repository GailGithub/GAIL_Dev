%% 1 Function and 2 CVs 
    w.func = @(x)[10*x(:,1)-5*x(:,2).^2+1*x(:,3).^3, x(:,1), x(:,2).^2];
    w.cv = [8,32/3]; hyperbox= [zeros(1,3);2*ones(1,3)];
    q = cubLattice_E(w,hyperbox,'uniform',1e-6,0); 

    %% 2 Functions and 0 CVs
    f = @(x)[sin(x),cos(x)-sin(1)+1-cos(1)]; hyperbox = [0;1];
    q = cubLattice_E(f,hyperbox,'uniform',1e-3,1e-2,'transform','C1');

    %% 2 Functions and 2 CVs
    w.func=@(x)[sin(x),cos(x)-sin(1)+1-cos(1) ,x-1/2,x.^2-1/3]
    w.cv = [0,0]; hyperbox= [0;1];
    q = cubSobol_g(w,hyperbox,'uniform',1e-6,0);

    %% 1 Function and 2 CVs 
    w.func=@(x)[sin(x), x-1/2,x.^2-1/3]
    w.cv = [0,0]; hyperbox= [0;1];
    q = cubSobol_g(w,hyperbox,'uniform',1e-6,0);

    %% 2 Functions and 3 CVs
    w.func=@(x)[sin(x),cos(x)-sin(1)+1-cos(1) ,x , x.^2, x.^3]
    w.cv = [0.5, 1/3, 1/4]; hyperbox= [0;1];
    q = cubLattice_E(w,hyperbox,'uniform',1e-6,0);

    %% 3 Functions and 2 CVs
    w.func=@(x)[sin(x),cos(x)-sin(1)+1-cos(1), sin(2.*x)/2.*cos(x), x-1/2,x.^2-1/3]
    w.cv = [0,0]; hyperbox= [0;1];
    q = cubLattice_E(w,hyperbox,'uniform',1e-6,0);

    %% 2 Functions and 1 CV 
    w.func=@(x)[x(:,1).^2+x(:,2).^2,x(:,1)+x(:,2)-1/3,(x(:,1))+x(:,2).^2]
    w.cv=[5/6 ];hyperbox=[0 0; 1 1]
    q = cubLattice_Kan(w,hyperbox,'uniform',1e-6,0, 2, 1);

    w.func=@(x)[x(:,1).^2+x(:,2).^2,x(:,1)+x(:,2)-1/3,(x(:,1))+x(:,2).^2]
    w.cv=[5/6];hyperbox=[0 0; 1 1]
    z = cubLattice_E(w,hyperbox,'uniform',1e-6,0, 2, 1);

%% Multi-dimensional (Using Keister's equation) 
gail.InitializeWorkspaceDisplay %initialize the workspace and the display parameters
normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
f1 = @(normt,a,d) ((2*pi*a^2).^(d/2)) * cos(a*sqrt(normt)) ...
   .* exp((1/2-a^2)*normt);
f = @(t,a,d) f1(normsqd(t),a,d);

abstol = 0; %absolute error tolerance
reltol = 0.01; %relative error tolerance
dvec = 1:5; %vector of dimensions
avec = [1 1/sqrt(2) 1/sqrt(1.5) 1/sqrt(3)]; %default value of a 
IMCvec = zeros(size(dvec)); %vector of answers
IMCvecEtime = zeros(size(dvec)); %vector of answers
IMCvecE = zeros(size(dvec)); %vector of answers

totE=[];
totTime=[];
for i=1:100
    for d = dvec
        w=@(t)[f(t, avec(1),d),f(t,avec(2),d),f(t,avec(3),d),f(t,avec(4),d)]
        hyperbox=[zeros(1,d); ones(1,d)];
        [IMCvecE(d), out, IMCvecEtime(d), ~] = cubLattice_E(w,hyperbox,'normal',1e-6,0);
    end 
    totE=vertcat(totE, IMCvecE);
    totTime=vertcat(totTime, IMCvecEtime);
end 

totE
totTime

%% 
totMCL=[];
totMCLTime=[];

for i=100 
    for d=dvec 
        
    end 
    
end 




 
%% 2 function, 1 CV
% dim 1
d=1
w.func=@(f(t, avec(1),d),f(t,avec(2),d), ]
w.cv = [0.5, 1/3, 1/4]; hyperbox=[zeros(1,d); ones(1,d)];

% dim 2
d=2
w.func=@(f(t, avec(1),d),f(t,avec(2),d), ]
w.cv = [0.5, 1/3, 1/4]; hyperbox=[zeros(1,d); ones(1,d)];

% dim 3 
d=2
w.func=@(f(t, avec(1),d),f(t,avec(2),d), ]
w.cv = [0.5, 1/3, 1/4]; hyperbox=[zeros(1,d); ones(1,d)];








