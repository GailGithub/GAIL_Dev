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

runs=1
totE=zeros(runs,5);
totTime=zeros(runs,5);
for i=1:runs
    for d = dvec
        if d==1;
            w.func=@(t)[f(t,1,1),f(t,1/sqrt(2),1),cos(t)]; %0.01
            w.cv=[1/sqrt(exp(1))];
            hyperbox=[zeros(1,d); ones(1,d)];
            [IMCvecE(d), out, IMCvecEtime(d), ~, ~] = cubLattice_E(w,hyperbox,'normal',1e-6,0);    
            
        elseif d==2;
            w.func=@(t)[f(t,1,2),f(t,1/sqrt(2),2), cos(t(:,1)).*cos(t(:,2))]; %0.02
            w.cv=[0.3677];
            hyperbox=[zeros(1,d); ones(1,d)];
            [IMCvecE(d), out, IMCvecEtime(d), ~] = cubLattice_E(w,hyperbox,'normal',1e-6,0);
            
        elseif d==3;
            w.func=@(t)[f(t,1,3),f(t,1/sqrt(2),3),cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3))]; % 0.03
            w.cv=[0.2231];
            hyperbox=[zeros(1,d); ones(1,d)];
            [IMCvecE(d), out, IMCvecEtime(d), ~] = cubLattice_E(w,hyperbox,'normal',1e-6,0);
            
        elseif d==4;
            w.func=@(t)[f(t,1,4),f(t,1/sqrt(2),4),cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3)).*cos(t(:,4))];%0.06
            w.cv=[0.1354];
            hyperbox=[zeros(1,d); ones(1,d)];
            [IMCvecE(d), out, IMCvecEtime(d), ~] = cubLattice_E(w,hyperbox,'normal',1e-6,0);
            
        elseif d==5;
            w.func=@(t)[f(t,1,5),f(t,1/sqrt(2),5),cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3)).*cos(t(:,4)).*cos(t(:,5))];  %0.5
            w.cv=[0.0821];
            hyperbox=[zeros(1,d); ones(1,d)];
            [IMCvecE(d), out, IMCvecEtime(d), ~] = cubLattice_E(w,hyperbox,'normal',1e-6,0);
        end
        
    end
    
    totE(i,:)=IMCvecE;
    totTime(i,:)=IMCvecEtime;
end

totE
totTime

