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
        if d==1;
            w.func=@(t)[f(t,1,1),f(t,1/sqrt(2),1),cos(t)]; %0.01
            YXn=@(n)f2(randn(n,1));
            a=@(t)cos(t)
            b=@(n)a(randn(n,1))
            c=meanMC_CLTKATE(b,0.001)
            w.cv=[c];
            hyperbox=[zeros(1,d); ones(1,d)];
            [IMCvecE(d), out, IMCvecEtime(d), ~] = cubLattice_E(w,hyperbox,'normal',1e-6,0);
            
        elseif d==2;
            f2=@(t)[f(t,1,2),f(t,1/sqrt(2),2),cos(t(:,1)).*cos(t(:,2))]; %0.02
            YXn=@(n)f2(randn(n,1));
            a=@(t)cos(t)
            b=@(n)a(randn(n,1))
            c=meanMC_CLTKATE(b,0.001)
            w.cv=[c];
            hyperbox=[zeros(1,d); ones(1,d)];
            [IMCvecE(d), out, IMCvecEtime(d), ~] = cubLattice_E(w,hyperbox,'normal',1e-6,0);
            
        elseif d==3;
            f2=@(t)[f(t,1,3),f(t,1/sqrt(2),3),cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3))]; % 0.03
            YXn=@(n)f2(randn(n,1));
            a=@(t)cos(t)
            b=@(n)a(randn(n,1))
            c=meanMC_CLTKATE(b,0.001)
            w.cv=[c];
            hyperbox=[zeros(1,d); ones(1,d)];
            [IMCvecE(d), out, IMCvecEtime(d), ~] = cubLattice_E(w,hyperbox,'normal',1e-6,0);
            
        elseif d==4;
            f2=@(t)[f(t,1,4),f(t,1/sqrt(2),4),cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3)).*cos(t(:,4))];%0.06
            YXn=@(n)f2(randn(n,1));
            a=@(t)cos(t)
            b=@(n)a(randn(n,1))
            c=meanMC_CLTKATE(b,0.001)
            w.cv=[c];
            hyperbox=[zeros(1,d); ones(1,d)];
            [IMCvecE(d), out, IMCvecEtime(d), ~] = cubLattice_E(w,hyperbox,'normal',1e-6,0);
            
        elseif d==5;
            
            f2=@(t)[f(t,1,5),f(t,1/sqrt(2),5),cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3)).*cos(t(:,4)).*cos(t(:,5))];  %0.5
            YXn=@(n)f2(randn(n,1));
            a=@(t)cos(t)
            b=@(n)a(randn(n,1))
            c=meanMC_CLTKATE(b,0.001)
            w.cv=[c];
            hyperbox=[zeros(1,d); ones(1,d)];
            [IMCvecE(d), out, IMCvecEtime(d), ~] = cubLattice_E(w,hyperbox,'normal',1e-6,0);
        end
        
    end
    totE=vertcat(totE, IMCvecE);
    totTime=vertcat(totTime, IMCvecEtime);
end

totE
totTime










