%% Example 7 doc test 
w.func = @(x)[10*x(:,1)-5*x(:,2).^2+1*x(:,3).^3, x(:,1), x(:,2).^2];
w.cv = [8,32/3]; hyperbox= [zeros(1,3);2*ones(1,3)];
q = cubLattice_E(w,hyperbox,'uniform',1e-6,0,1,2); exactsol = 128/3; 


%% 1 Function and 2 CVs 
w.func = @(x)[10*x(:,1)-5*x(:,2).^2+1*x(:,3).^3, x(:,1), x(:,2).^2];
w.cv = [8,32/3]; hyperbox= [zeros(1,3);2*ones(1,3)];
q = cubLattice_Kan(w,hyperbox,'uniform',1e-6,0,1,2); 

%% 2 Functions and 0 CVs
f = @(x)[sin(x),cos(x)-sin(1)+1-cos(1)]; hyperbox = [0;1];
q = cubLattice_Kan(f,hyperbox,'uniform',1e-3,1e-2,'transform','C1');

%% 2 Functions and 2 CVs
w.func=@(x)[sin(x),cos(x)-sin(1)+1-cos(1) ,x-1/2,x.^2-1/3]
w.cv = [0,0]; hyperbox= [0;1];
q = cubLattice_E(w,hyperbox,'uniform',1e-6,0, 2, 2);

%% 1 Function and 2 CVs 
w.func=@(x)[sin(x), x-1/2,x.^2-1/3]
w.cv = [0,0]; hyperbox= [0;1];
q = cubLattice_Kan(w,hyperbox,'uniform',1e-6,0, 1, 2);

%% 2 Functions and 3 CVs
w.func=@(x)[sin(x),cos(x)-sin(1)+1-cos(1) ,x , x.^2, x.^3]
w.cv = [0.5, 1/3, 1/4]; hyperbox= [0;1];
q = cubLattice_E(w,hyperbox,'uniform',1e-6,0, 2, 3);

%% 3 Functions and 2 CVs
w.func=@(x)[sin(x),cos(x)-sin(1)+1-cos(1), sin(2.*x)/2.*cos(x), x-1/2,x.^2-1/3]
w.cv = [0,0]; hyperbox= [0;1];
q = cubLattice_E(w,hyperbox,'uniform',1e-6,0, 3, 2);






