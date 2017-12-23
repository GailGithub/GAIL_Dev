%Run TestcubMC on the step function
function [res,test,fun,param]=RunTestcubMConGaussiand1()
format compact

%test.nrep=500;%in the paper, we use 500 repilcation numbers
test.nrep = 500;
fun.funtype='gaussian';
param.dim=1;
param.measure='uniform';

% param.measure='normal';
% param.interval=[-inf(1,param.dim);inf(1,param.dim)];
param.impyes=false;
param.abstol=1e-3;
param.n0=2^13;
param.interval=[zeros(1,param.dim);ones(1,param.dim)];
param.bmina=param.interval(2,:)-param.interval(1,:);
test.howoftenrep=10;
shapemin=1e-6;
shapemax=1;

test.randch.shapeoverall=...
    shapemin*(shapemax/shapemin).^rand(test.nrep,param.dim);
scalemin=.1;
scalemax=10;

test.randch.scaleoverall=...
    scalemin*(scalemax/scalemin).^rand(test.nrep,param.dim);
test.randch.centercoverall=rand(test.nrep,param.dim);

overmultcmin=1e-2;
overmultcmax=1e2; 
test.randch.addcoverall=ones(test.nrep,param.dim);

test.randch.gaussian1=sqrt(pi)*test.randch.shapeoverall./(2*repmat(param.bmina,test.nrep,1)).*...
        (erf((repmat(param.interval(2,:),test.nrep,1)-test.randch.centercoverall)...
        ./test.randch.shapeoverall)-erf((repmat(param.interval(1,:),test.nrep,1)-...
        test.randch.centercoverall)./test.randch.shapeoverall));
test.randch.gaussian2=sqrt(pi)*test.randch.shapeoverall./...
    (2*sqrt(2)*repmat(param.bmina,test.nrep,1))...
        .*(erf(sqrt(2)*(repmat(param.interval(2,:),test.nrep,1) ...
        -test.randch.centercoverall)...
        ./test.randch.shapeoverall) - ...
        erf(sqrt(2)*(repmat(param.interval(1,:),test.nrep,1) ...
        -test.randch.centercoverall)...
        ./test.randch.shapeoverall));    
test.randch.overmultcoverall=sqrt(overmultcmin*(overmultcmax/overmultcmin)...
    .^rand(test.nrep,1)./(prod((test.randch.addcoverall.^2+...
    2*test.randch.addcoverall.*test.randch.scaleoverall.*test.randch.gaussian1...
    +test.randch.scaleoverall.^2.*test.randch.gaussian2),2)-...
    prod((test.randch.addcoverall+test.randch.scaleoverall.* ...
    test.randch.gaussian1).^2,2)));
test.randch.overaddcoverall=ones(test.nrep,1)-...
    (test.randch.overmultcoverall.*prod((test.randch.addcoverall ...
        + test.randch.scaleoverall.*test.randch.gaussian1),2));


test.randchoicefun=@randchoicegaussiand1;
test.whichsample={'iid','iidheavy'};
if exist('chebfun','class')==8 && exist('sobolset','class')==8
    test.whichsample={'iid','iidheavy','Sobol','Sobolheavy','quad','quadgk','chebfun','chebfunheavy'};
elseif exist('chebfun','class')==8 && exist('sobolset','class')~= 8
    test.whichsample={'iid','iidheavy','quad','quadgk','chebfun','chebfunheavy'};
elseif exist('chebfun','class')~=8 && exist('sobolset','class')==8
    test.whichsample={'iid','iidheavy','Sobol','Sobolheavy','quad','quadgk',};
else
    test.whichsample={'iid','iidheavy','quad','quadgk'};
end
%test.whichsample={'iid','iidheavy','Sobol','Sobolheavy'};
res = TestcubMCDiffSettings(test,fun,param);
end
