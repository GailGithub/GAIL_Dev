%a = 1:5;
normsqd = @(t) sum(t.*t,2);
f1 = @(normt,a,d) ((2*pi*a^2).^(d/2)) * cos(a*sqrt(normt)) ...
   .* exp((1/2-a^2)*normt);
f = @(t,a,d) f1(normsqd(t),a,d);
g=100
Val=zeros(1,g);
Time=zeros(1,g);
for i=1:g
%f2=@(t)[f(t,1,1),f(t,1/sqrt(2),1),cos(t)];
%dimension1
f2=@(t)[f(t,1,1),f(t,1/sqrt(2),1)];% 0.01
f2=@(t)[f(t,1,1),f(t,1/sqrt(2),1),cos(t)]; %0.01
YXn=@(n)f2(randn(n,1));
a=@(t)cos(t)
b=@(n)a(randn(n,1))
c=meanMC_CLTKATE(b,0.001)
s=struct('Y',YXn,'nY',2,'trueMuCV',c)
[Val(i), ~, Time]= meanMC_CLTKATE(s,0,0.01);

end





% %dimension2
% f2=@(t)[f(t,1,2),f(t,1/sqrt(2),2)];%creat a matrix with same mean for each column  0.03
% f2=@(t)[f(t,1,2),f(t,1/sqrt(2),2),cos(t(:,1)).*cos(t(:,2))]; %0.02
% YXn=@(n)f2(randn(n,2));
% a=@(t)cos(t(:,1)).*cos(t(:,2))
% b=@(n)a(randn(n,2))
% c=meanMC_CLTKATE(b,0.001)
% s=struct('Y',YXn,'nY',2,'trueMuCV',c)
% [hmu, out]= meanMC_CLTKATE(s,0,0.01);
% hmu
% out

%dimension3
% f2=@(t)[f(t,1,3),f(t,1/sqrt(2),3)];%creat a matrix with same mean for each column 0.05
% f2=@(t)[f(t,1,3),f(t,1/sqrt(2),3),cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3))]; % 0.03
% YXn=@(n)f2(randn(n,3));
% a=@(t)cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3))
% b=@(n)a(randn(n,3))
% c=meanMC_CLTKATE(b,0.001)
% s=struct('Y',YXn,'nY',2,'trueMuCV',c)
% [hmu, out]= meanMC_CLTKATE(s,0,0.01);
% hmu
% out

%dimension4
% f2=@(t)[f(t,1,4),f(t,1/sqrt(2),4)];%creat a matrix with same mean for each column 0.12
% f2=@(t)[f(t,1,4),f(t,1/sqrt(2),4),cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3)).*cos(t(:,4))];%0.06
% YXn=@(n)f2(randn(n,4));
% a=@(t)cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3)).*cos(t(:,4))
% b=@(n)a(randn(n,4))
% c=meanMC_CLTKATE(b,0.001)
% s=struct('Y',YXn,'nY',2,'trueMuCV',c)
% [hmu, out]= meanMC_CLTKATE(s,0,0.01);
% hmu
% out

%dimension5
% f2=@(t)[f(t,1,5),f(t,1/sqrt(2),5)];%creat a matrix with same mean for each column 1
% f2=@(t)[f(t,1,5),f(t,1/sqrt(2),5),cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3)).*cos(t(:,4)).*cos(t(:,5))];  %0.5
% YXn=@(n)f2(randn(n,5));
% a=@(t)cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3)).*cos(t(:,4)).*cos(t(:,5))
% b=@(n)a(randn(n,5))
% c=meanMC_CLTKATE(b,0.001)
% s=struct('Y',YXn,'nY',2,'trueMuCV',c)
% [hmu, out]= meanMC_CLTKATE(s,0,0.01);
% out
% hmu





% %f2=@(t)[f(t,1,3),f(t,1/sqrt(2),3),cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3))];
% YXn=@(n)f2(randn(n,3));
% %abstol=0;
% %reltol=1e-3;
% tic
% %s=struct('Y',YXn,'nY',2,'trueMuCV',[sin(1).^3 3.*sin(1)]); %create a structure as input for meanMC_CLT
% %s=struct('Y',YXn,'nY',2,'trueMuCV',sin(1)*sin(1)+sin(1))  0.2219 0.6064
% %0.3674
% s=struct('Y',YXn,'nY',2)
% [hmu, out]= meanMC_CLTKATE(s,0,0.01);
% toc

%%
YX=@(x)[sin(x),cos(x)-sin(1)+1-cos(1),x-1/2,x.^2-1/3]
%YX=@(x)[sin(x),cos(x)-sin(1)+1-cos(1)]
YXn=@(n)YX(rand(n,1))
s=struct('Y',YXn,'nY',2,'trueMuCV',zeros(2,1)')
%s=struct('Y',YXn,'nY',2)
[hmu,mean_out]=meanMC_CLTKATE(s,0.01)

%%
f = @(x) [x(:,1).^3.*x(:,2).^3.*x(:,3).^3, x(:,1).*x(:,2).*x(:,3)];
%f = @(x)x(:,1)+x(:,2)+x(:,3);
Yrand=@(n)f(rand(n,3));
s=struct('Y',@(n)f(rand(n,3)),'nY',1)
[hmu,mean_out]=meanMC_CLTKATE(Yrand,1e-4)


%%
f = @(x) [x(:,1).^3.*x(:,2).^3.*x(:,3).^3, x(:,1).^2.*x(:,2).^2.*x(:,3).^2-1/27+1/64,x(:,1).*x(:,2).*x(:,3),x(:,1)+x(:,2)+x(:,3)];
%f = @(x) [x(:,1).^3.*x(:,2).^3.*x(:,3).^3];
s=struct('Y',@(n)f(rand(n,3)),'nY',2,'trueMuCV',[1/8 1.5])
[hmu,mean_out]=meanMC_CLTKATE(s,1e-4)

%%
%f = @(x) [x(:,1).^3.*x(:,2).^3.*x(:,3).^3, x(:,1).^2.*x(:,2).^2.*x(:,3).^2-1/27+1/64,x(:,1).*x(:,2).*x(:,3),x(:,1)+x(:,2)+x(:,3)];
f = @(x) (cos(x(:,1))+cos(x(:,2))+cos(x(:,3)));
Yrand=@(n)f(rand(n,3));
%s=struct('Y',@(n)f(rand(n,3)),'nY',2,'trueMuCV',[1/8 1.5])
[hmu,mean_out]=meanMC_CLTKATE(Yrand,1e-4)

%%
f1 = @(normt,a,d) ((2*pi*a^2).^(d/2)) * cos(a*sqrt(normt)) ...
   .* exp((1/2-a^2)*normt);
f = @(t,a,d) f1(normsqd(t),a,d);
abstol = 0; %absolute error tolerance
reltol = 0.1; %relative error tolerance
%dvec = 1:5; %vector of dimensions
avec = [0.5 1/sqrt(2)]; %default value of a 



return          
            
            


obj2=obj;
            obj2.payoffParam.optType={'gmean'};
            YX1n=@(n)genOptPayoffs(obj,n);
            YX2n=@(n)genOptPayoffs(obj2,n);
            YX=[YX1n, YX2n];
            s=struct('Yrand',YX,'q',1,'xMean',2.168443483018763);
            [price, outtemp] = meanMC_CLT(s, ...
               obj.priceParam.absTol, obj.priceParam.relTol, ...
               obj.priceParam.alpha);
            
            
            
%% Example 1:
% Estimate the integral with integrand f(x) = x1.*x2 in the interval [0,1)^2 with absolute 
% tolerance 1e-5 and relative tolerence 0:
% 
 f = @(x) prod(x,2);
 q = meanMC_CLT(@(n)f(rand(n,2)),1e-5,0); exactsol = 1/4; 
 check = abs(exactsol-q) < 1e-5
 
 %%
%check = 1
%
%
% Example 2:
% Estimate the integral with integrand f(x) = x1.^3.*x2.^3.*x3.^3
% in the interval [0,1)^3 with pure absolute error 1e-5 using x1.*x2.*x3 as control variate:
% 
f=@(x) [x(:,1).^3.*x(:,2).^3.*x(:,3).^3, x(:,1).*x(:,2).*x(:,3)];
 s=struct('Y',@(n)f(rand(n,3)),'nY',1,'trueMuCV',1/8)
 [hmu,mean_out]=meanMC_CLTKATE(s,1e-3,0); exactsol = 1/64;
 check = abs(exactsol-hmu) < max(1e-3,1e-3*abs(exactsol))
 
 %%
% check = 1
%
% Example 3:
% Estimate the integrals with integrands f1(x) = x1.^3.*x2.^3.*x3.^3 and 
% f2(x)= x1.^2.*x2.^2.*x3.^2-1/27+1/64 in the interval [0,1)^3
% using  x1.*x2.*x3 and x1+x2.^3+x3 as control variate:
 f = @(x) [x(:,1).^3.*x(:,2).^3.*x(:,3).^3, x(:,1).^2.*x(:,2).^2.*x(:,3).^2-1/27+1/64,x(:,1).*x(:,2).*x(:,3),x(:,1)+x(:,2)+x(:,3)];
 s=struct('Y',@(n)f(rand(n,3)),'nY',2,'trueMuCV',[1/8 1.5])
 [hmu,mean_out]=meanMC_CLTKATE(s,1e-4,1e-3); exactsol = 1/64;
 check = abs(exactsol-hmu) < max(1e-4,1e-3*abs(exactsol))
% check = 1

%%
normsqd = @(t) sum(t.*t,2);
f1 = @(normt,a,d) ((2*pi*a^2).^(d/2)) * cos(a*sqrt(normt)) ...
   .* exp((1/2-a^2)*normt);
f = @(t,a,d) f1(normsqd(t),a,d);
dim=2
f2=@(t)[f(t,1,dim),f(t,1/sqrt(2),dim),f(t,1/sqrt(1.5),dim),f(t,1/sqrt(3),dim)]
YXn=@(n)f2(randn(n,1));            
abstol=0;
reltol=1e-3;
tic
%s=struct('Y',YXn,'nY',2,'trueMuCV',[sin(1).^3 3.*sin(1)]); %create a structure as input for meanMC_CLT
%s=struct('Y',YXn,'nY',2,'trueMuCV',sin(1)*sin(1)+sin(1)) 
s=struct('Y',YXn,'nY',2,'trueMuCV',[1.8045 1.8045])
[hmu, out]= meanMC_CLTKATE(s,0,0.01);
toc          


%%

abstol = 0; %absolute error tolerance
reltol = 0.01; %relative error tolerance
dvec = 1:5; %vector of dimensions
avec = [1 1/sqrt(2) 1/sqrt(1.5) 1/sqrt(3)]; %default value of a 
IMCvec = zeros(size(dvec)); %vector of answers
f2= @(t,d) cell2mat(arrayfun(@(a) f(t,a,d),avec,'UniformOutput',false));

 for d = dvec
 f3=@(t)f2(t,d);
 YXn=@(n)f3(randn(n,d));
 s=struct('Y',YXn,'nY',size(avec,2)); %create a structure as input for meanMC_CLT
 tic
 [IMCvec(d),out]= meanMC_CLTKATE(s,abstol,reltol);
 toc
 display(out)
 end
%%
normsqd = @(t) sum(t.*t,2);
f1 = @(normt,a,d) ((2*pi*a^2).^(d/2)) * cos(a*sqrt(normt)) ...
   .* exp((1/2-a^2)*normt);
f = @(t,a,d) f1(normsqd(t),a,d);
% %f2=@(t)[f(t,1,1),f(t,1/sqrt(2),1),cos(t)];
 
%% Dimension 1 
g=100
Val1=zeros(1,g);
Time1=zeros(1,g);
for i=1:g
    %f2=@(t)[f(t,1,1),f(t,1/sqrt(2),1)];% 0.01
     f2=@(t)[f(t,1,1),f(t,1/sqrt(2),1),cos(t)]; %0.01
    YXn=@(n)f2(randn(n,1));
    a=@(t)cos(t)
    b=@(n)a(randn(n,1))
    c=meanMC_CLTKATE(b,0.0001)
    s=struct('Y',YXn,'nY',2,'trueMuCV',c)
    [Val1(i), ~, Time1(i)]= meanMC_CLTKATE(s,0,0.01);
end 

%%
%%Example 4
normsqd = @(t) sum(t.*t,2);
f = @(normt,a,d) ((2*pi*a^2).^(d/2)) * cos(a*sqrt(normt)) ...
   .* exp((1/2-a^2)*normt);
f1 = @(t,a,d) f(normsqd(t),a,d);
     f2=@(t)[f1(t,1,1),f1(t,1/sqrt(2),1),cos(t)]; %0.01
    YXn=@(n)f2(randn(n,1));
    s=struct('Y',YXn,'nY',2,'trueMuCV',1/sqrt(exp(1)))
    [hmu,out]= meanMC_CLTKATE(s,0,0.0001);
    
%% Dimension2
 
g=100
Val2=zeros(1,g);
Time2=zeros(1,g);
for i=1:g
    %f2=@(t)[f(t,1,2),f(t,1/sqrt(2),2)];%creat a matrix with same mean for each column  0.03
     f2=@(t)[f(t,1,2),f(t,1/sqrt(2),2),cos(t(:,1)).*cos(t(:,2))]; %0.02
    YXn=@(n)f2(randn(n,2));
    a=@(t)cos(t(:,1)).*cos(t(:,2))
    b=@(n)a(randn(n,2))
    c=meanMC_CLTKATE(b,0.001)
    s=struct('Y',YXn,'nY',2,'trueMuCV',c)
    [Val2(i), ~, Time2(i)]= meanMC_CLTKATE(s,0,0.01);
end 
 



%% Dimension3
g=100
Val3=zeros(1,g);
Time3=zeros(1,g);
for i=1:g
    %f2=@(t)[f(t,1,3),f(t,1/sqrt(2),3)];%creat a matrix with same mean for each column 0.05
     f2=@(t)[f(t,1,3),f(t,1/sqrt(2),3),cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3))]; % 0.03
    YXn=@(n)f2(randn(n,3));
    a=@(t)cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3))
    b=@(n)a(randn(n,3))
    c=meanMC_CLTKATE(b,0.001)
    s=struct('Y',YXn,'nY',2,'trueMuCV',c)
    [Val3(i), ~, Time3(i)]= meanMC_CLTKATE(s,0,0.01);
end
            

%% Dimension 4
g=100
Val4=zeros(1,g);
Time4=zeros(1,g);
for i=1:g
    f2=@(t)[f(t,1,4),f(t,1/sqrt(2),4)];%creat a matrix with same mean for each column 0.12
    % f2=@(t)[f(t,1,4),f(t,1/sqrt(2),4),cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3)).*cos(t(:,4))];%0.06
    YXn=@(n)f2(randn(n,4));
    a=@(t)cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3)).*cos(t(:,4))
    b=@(n)a(randn(n,4))
    c=meanMC_CLTKATE(b,0.001)
    s=struct('Y',YXn,'nY',2,'trueMuCV',c)
    [Val4(i), ~, Time4(i)]= meanMC_CLTKATE(s,0,0.01);
end


%% Dimension 5
g=100
Val5=zeros(1,g);
Time5=zeros(1,g);
for i=1:g
    f2=@(t)[f(t,1,5),f(t,1/sqrt(2),5)];%creat a matrix with same mean for each column 1
    %f2=@(t)[f(t,1,5),f(t,1/sqrt(2),5),cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3)).*cos(t(:,4)).*cos(t(:,5))];  %0.5
    YXn=@(n)f2(randn(n,5));
    a=@(t)cos(t(:,1)).*cos(t(:,2)).*cos(t(:,3)).*cos(t(:,4)).*cos(t(:,5))
    b=@(n)a(randn(n,5))
    c=meanMC_CLTKATE(b,0.0001)
    s=struct('Y',YXn,'nY',2,'trueMuCV',c)
    [Val5(i), ~,Time5(i)]= meanMC_CLTKATE(s,0,0.01);
end