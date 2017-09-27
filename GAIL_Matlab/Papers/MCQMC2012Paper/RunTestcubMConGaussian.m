%Run TestcubMC on the step function
function [res,test,fun,param] = RunTestcubMConGaussian()
format compact
%test.nrep=500; in the paper, we use 500 repilcation numbers
test.nrep = 100;
fun.funtype='gaussian';
param.dim=8;
param.measure='uniform';

% param.measure='normal';
% param.interval=[-inf(1,param.dim);inf(1,param.dim)];
dimchoice=[2 3 4 5 6 7 8]';
ndim=size(dimchoice,1);
test.randch.dimoverall=dimchoice(randi(ndim,test.nrep,1));

param.impyes=false;
param.abstol=1e-2;
param.n0=2^13;
test.howoftenrep=10;
shapemin=1e-6;
shapemax=1;
test.randch.shapeoverall=nan(test.nrep,max(dimchoice));
scalemin=.1;
scalemax=10;
test.randch.scaleoverall=nan(test.nrep,max(dimchoice));
test.randch.centercoverall=nan(test.nrep,max(dimchoice));
overmultcmin=1e-2;
overmultcmax=1e2; 
test.randch.addcoverall=nan(test.nrep,max(dimchoice));
param.interval=nan(2*test.nrep,max(dimchoice));
param.bmina=nan(test.nrep,max(dimchoice));
param.interval1=nan(test.nrep,max(dimchoice));
param.interval2=nan(test.nrep,max(dimchoice));
test.randch.overaddcoverall=nan(test.nrep,1);
test.randch.overmultcoverall=nan(test.nrep,1);
for i=1:test.nrep
    param.interval(2*i-1:2*i,1:test.randch.dimoverall(i))= ...
        [zeros(1,test.randch.dimoverall(i));ones(1,test.randch.dimoverall(i))];
    
    param.interval1(i,1:test.randch.dimoverall(i))= ...
        zeros(1,test.randch.dimoverall(i));
    
    param.interval2(i,1:test.randch.dimoverall(i))= ...
        ones(1,test.randch.dimoverall(i));
    
    param.bmina(i,1:test.randch.dimoverall(i))=...
        param.interval(2*i,test.randch.dimoverall(i))-...
        param.interval(2*i-1,test.randch.dimoverall(i));
    
    test.randch.shapeoverall(i,1:test.randch.dimoverall(i))...
        =shapemin*(shapemax/shapemin).^rand(1,test.randch.dimoverall(i));

    test.randch.scaleoverall(i,1:test.randch.dimoverall(i))...
        =scalemin*(scalemax/scalemin).^rand(1,test.randch.dimoverall(i));

    test.randch.centercoverall(i,1:test.randch.dimoverall(i))...
        =rand(1,test.randch.dimoverall(i));

    test.randch.addcoverall(i,1:test.randch.dimoverall(i))...
        =ones(1,test.randch.dimoverall(i));
end

test.randch.gaussian1=sqrt(pi)*test.randch.shapeoverall./...
    (2*param.bmina).*...
        (erf((param.interval2-test.randch.centercoverall)...
        ./test.randch.shapeoverall)-erf((param.interval1-...
        test.randch.centercoverall)./test.randch.shapeoverall));
test.randch.gaussian2=sqrt(pi)*test.randch.shapeoverall./...
    (2*sqrt(2)*param.bmina)...
        .*(erf(sqrt(2)*(param.interval2 ...
        -test.randch.centercoverall)...
        ./test.randch.shapeoverall) - ...
        erf(sqrt(2)*(param.interval1 ...
        -test.randch.centercoverall)...
        ./test.randch.shapeoverall));
    test.randch.moment1=test.randch.addcoverall+test.randch.scaleoverall.* ...
    test.randch.gaussian1;
    test.randch.moment1sq=(test.randch.addcoverall+test.randch.scaleoverall.* ...
    test.randch.gaussian1).^2;
test.randch.moment2=(test.randch.addcoverall.^2+...
    2*test.randch.addcoverall.*test.randch.scaleoverall.*test.randch.gaussian1...
    +test.randch.scaleoverall.^2.*test.randch.gaussian2);
for i=1:test.nrep
test.randch.overmultcoverall(i)=sqrt(overmultcmin*(overmultcmax/overmultcmin)...
    .^rand(1)./(prod(test.randch.moment2(i,1:test.randch.dimoverall(i)))-...
    prod(test.randch.moment1sq(i,1:test.randch.dimoverall(i)))));
test.randch.overaddcoverall(i)=1-...
    (test.randch.overmultcoverall(i).*prod(test.randch.moment1(i,1:test.randch.dimoverall(i))));
end

test.randchoicefun=@randchoicegaussian;

%test.whichsample={'iid','iidheavy','Sobol','Sobolheavy','quad','quadgk','chebfun','chebfunheavy'};
test.whichsample={'iid','iidheavy','Sobol','Sobolheavy'};
res = TestcubMCDiffSettings(test,fun,param);
end
