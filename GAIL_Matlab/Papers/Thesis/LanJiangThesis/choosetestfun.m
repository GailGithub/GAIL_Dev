function [testfun,param]=choosetestfun(fun,param)
%   This function chooses and sets up a test function from the parameters
%      input by the user and contained in the structures
%      fun and param
%   fun.funtype         = type of test function
%   param.interval      = domain of test function
%   param.dim           = dimension of the domain
%   param.measure           = probability density function for integration
%   param.exactintegral = exact value of the integral (scalar)
%   fun.shape           = shape parameter (1 x param.dim)
%   fun.scale           = scale parameter (1 x param.dim)
%   fun.addc            = additive constant (1 x param.dim)
%   fun.overaddc        = overall additive constant (scalar)
%   fun.overmultc       = overall multiplicative constant (scalar)

if nargin < 2; %give the basic default parameters 
    param.interval=[0;1]; %default integration interval
    if nargin < 1; fun.funtype='exp'; end %exponential test function
end
if ~isfield(param,'interval'); param.interval=[0;1]; end %default interval

[~,param]=cubMCparam([],param,'fun'); %check validity of some parameters

if strcmp(fun.funtype,'exp') %exponential test function
    [testfun,param]=makeExpTestFun(fun,param);
elseif strcmp(fun.funtype,'step') %square step test function
    [testfun,param]=makeStepTestFun(fun,param);
elseif strcmp(fun.funtype,'gaussian')
    [testfun,param]=makeGaussianTestFun(fun,param);
elseif strcmp(fun.funtype,'gaussianker')
     [testfun,param]=makeGaussianKerTestFun(fun,param);
else
    error('Function type not recognized')
end

param.fun=fun; %copy of the function parameters
end

%% Exponential Test Function
function [testfun,param]=makeExpTestFun(fun,param)
%Create the exponential test function
fun=verifyparam(fun,{'shape','scale','addc','overaddc','overmultc'}, ...
    {[1 param.dim],[1 param.dim],[1 param.dim],[1 1],[1 1]}, ...
    {1,1,0,0,1});
testfun=@(x) expfun(x,fun.overaddc,fun.overmultc,fun.addc,fun.scale,...
    fun.shape);

%% Compute exact integral of this function
if strcmp(param.measure,'uniform')
    bmina=param.interval(2,:)-param.interval(1,:);
    
    param.exactintegral=fun.overaddc*prod(bmina) + ...
        fun.overmultc.*prod(fun.addc.*bmina+fun.scale.*(exp(fun.shape ...
        .*param.interval(2,:))-exp(fun.shape.*param.interval(1,:)))./fun.shape);
    
    param.exactmean=fun.overaddc+ fun.overmultc.*prod(fun.addc ...
        + fun.scale.*(exp(fun.shape*diag(param.interval(2,:))) ...
        - exp(fun.shape*diag(param.interval(1,:))))./(fun.shape*diag(bmina)));
    
    param.moment1=fun.overmultc.*prod(fun.addc ...
        + fun.scale.*(exp(fun.shape*diag(param.interval(2,:))) ...
        - exp(fun.shape*diag(param.interval(1,:))))./(fun.shape*diag(bmina)));
    
    param.moment2=fun.overmultc.^2.*...
        prod(fun.addc.^2+2*fun.addc.*fun.scale.*...
        (exp(fun.shape*diag(param.interval(2,:))) ...
        - exp(fun.shape*diag(param.interval(1,:))))./(fun.shape*diag(bmina))+fun.scale.^2.*...
        (exp(2*fun.shape*diag(param.interval(2,:))) ...
        - exp(2*fun.shape*diag(param.interval(1,:))))./(2*fun.shape*diag(bmina)));
  
    param.exactvariance=param.moment2-param.moment1.^2;
    
    param.moment3=fun.overmultc.^3.*prod(fun.addc.^3+...
        3*fun.addc.^2.*fun.scale.*(exp(fun.shape*diag(param.interval(2,:))) ...
        - exp(fun.shape*diag(param.interval(1,:))))./(fun.shape*diag(bmina))+...
        3*fun.addc.*fun.scale.^2.*(exp(2*fun.shape*diag(param.interval(2,:))) ...
        - exp(2*fun.shape*diag(param.interval(1,:))))./(2*fun.shape*diag(bmina))+...
        fun.scale.^3.*(exp(3*fun.shape*diag(param.interval(2,:))) ...
        - exp(3*fun.shape*diag(param.interval(1,:))))./(3*fun.shape*diag(bmina)));
    
    param.moment4=fun.overmultc.^4.*prod(fun.addc.^4+...
        4*fun.addc.^3.*fun.scale.*(exp(fun.shape*diag(param.interval(2,:))) ...
        - exp(fun.shape*diag(param.interval(1,:))))./(fun.shape*diag(bmina))+...
        6*fun.addc.^2.*fun.scale.^2.*(exp(2*fun.shape*diag(param.interval(2,:))) ...
        - exp(2*fun.shape*diag(param.interval(1,:))))./(2*fun.shape*diag(bmina))+...
        4*fun.addc.*fun.scale.^3.*(exp(3*fun.shape*diag(param.interval(2,:))) ...
        - exp(3*fun.shape*diag(param.interval(1,:))))./(3*fun.shape*diag(bmina))+...
        fun.scale.^4.*(exp(4*fun.shape*diag(param.interval(2,:))) ...
        - exp(4*fun.shape*diag(param.interval(1,:))))./(4*fun.shape*diag(bmina)));
    
    param.exactkurtosis=(param.moment4-4*param.moment1...
        .*param.moment3+6*param.moment1.^2.*param.moment2...
        -4*param.moment1.^4+param.moment1.^4)...
        ./(param.exactvariance.^2);
    
    
else %normal distribution
    param.exactintegral=fun.overaddc ...
        + fun.overmultc.*prod(fun.addc ...
        + fun.scale.*exp(fun.shape.^2/2));
    exactintegralminus=fun.overmultc.*prod(fun.addc ...
        + fun.scale.*exp(fun.shape.^2/2));
    moment2=fun.overmultc.^2.*prod(fun.addc.^2 ...
        + 2.*fun.addc.*fun.scale.*exp(fun.shape.^2/2)...
        + fun.scale.^2.*exp((2*fun.shape).^2/2));
    moment3=fun.overmultc.^3.*prod(fun.addc.^3 ...
        + 3.*fun.addc.^2.*fun.scale.*exp(fun.shape.^2/2)...
        + 3*fun.addc.*fun.scale.^2.*exp((2*fun.shape).^2/2)...
        +fun.scale.^3.*exp((3*fun.shape).^2/2));
    moment4=fun.overmultc.^4.*prod(fun.addc.^4 ...
        + 4.*fun.addc.^3.*fun.scale.*exp(fun.shape.^2/2)...
        + 6*fun.addc.^2.*fun.scale.^2.*exp((2*fun.shape).^2/2)...
        +fun.addc.*fun.scale.^3.*exp((3*fun.shape).^2/2)...
        +fun.scale.^4.*exp(4*fun.shape).^2/2);
    param.exactvariance=moment2-exactintegralminus.^2;
      param.exactkurtosis=(moment4-4*exactintegralminus...
        .*moment3+6*exactintegralminus.^2.*moment2...
        -4*exactintegralminus.^4+exactintegralminus.^4)...
        ./(param.exactvariance.^2);
    
end

%% Prepare text description of the exponential function
if fun.overaddc==0; 
    formulastring='';
else
    formulastring=[num2str(fun.overaddc) ' + ']; 
end
if fun.overmultc~=1; 
    formulastring=[formulastring num2str(fun.overmultc) ' ']; 
end
if param.dim~=1; 
    formulastring=[formulastring 'prod_{j=1}^' int2str(param.dim)];
end
formulastring=[formulastring '['];
numdim=min(6,param.dim);
if param.dim<=numdim
    paramendstring=char(10); 
else
    paramendstring=[' ...' char(10)];
end
if fun.addc==0; 
    paramstring='';
else
    formulastring=[formulastring 'a_j + '];
    paramstring=['    a_j = ' num2str(fun.addc(1:numdim),'%6.3g')...
        paramendstring]; 
end
if fun.scale~=1; 
    formulastring=[formulastring 'b_j']; 
    paramstring=[paramstring ...
        '    b_j = ' num2str(fun.scale(1:numdim),'%6.3g')...
        paramendstring]; 
end
if fun.shape==1; 
    formulastring=[formulastring 'exp(x_j)'];
else
    formulastring=[formulastring 'exp(c_j x_j)'];
    paramstring=[paramstring ...
        '    c_j = ' num2str(fun.shape(1:numdim),'%6.3g')...
        paramendstring]; 
end
param.funDescribe=['   f(x) = ' formulastring ']' char(10) paramstring];
end

function f=expfun(x,overaddc,overmultc,addc,scale,shape)
    n=size(x,1);
    f=overaddc+overmultc.*prod(repmat(addc,n,1) ...
        + repmat(scale,n,1).*exp(repmat(shape,n,1).*x),2);
end

%% Step Test Function
function [testfun,param]=makeStepTestFun(fun,param)
%Create the exponential test function
fun=verifyparam(fun,{'shift','shape','scale','addc',...
    'overaddc','overmultc'}, ...
    {[1 param.dim],[1 param.dim],[1 param.dim],[1 param.dim],...
    [1 1],[1 1]}, ...
    {0.5,0.5,1,1,0,0});
bmina=param.interval(2,:)-param.interval(1,:);
testfun=@(x) stepfun(x,fun.overaddc,fun.overmultc,fun.addc,fun.scale,...
    fun.shape,fun.shift,bmina);

%% Compute exact integral of this function
%   as well as the variance and kurtosis
if strcmp(param.measure,'uniform')
    prodbmina=prod(bmina);
    moment1pc=fun.addc.*bmina+fun.scale.*fun.shape;
    moment2pc=(fun.addc.^2).*bmina+(2.*fun.addc+fun.scale)...
        .*fun.scale.*fun.shape;
    moment3pc=(fun.addc.^3).*bmina+(3.*(fun.addc.^2)...
        +fun.scale.*(3.*fun.addc+fun.scale)).*fun.scale.*fun.shape;
    moment4pc=(fun.addc.^4).*bmina+(4.*(fun.addc.^3)...
        +fun.scale.*(6.*(fun.addc.^2)+fun.scale.*(4.*fun.addc...
        +fun.scale))).*fun.scale.*fun.shape;
    
    moment1=prod(moment1pc);
    moment2=prod(moment2pc);
    moment3=prod(moment3pc);
    moment4=prod(moment4pc);
    param.exactintegral=fun.overaddc.*prodbmina + fun.overmultc.*moment1;
    param.exactvariance=(fun.overmultc.^2).*(moment2./prodbmina...
        -moment1.^2./(prodbmina)^2);
    param.exactkurtosis=(fun.overmultc.^4).*(...
        moment4./prodbmina+moment1.*(-4*moment3./(prodbmina)^2 ...
    +moment1.*(6*moment2./(prodbmina)^3 ...
        -3*moment1.^2./(prodbmina)^4)))./(param.exactvariance.^2);
    %keyboard
end

%% Prepare text description of the step function
if fun.overaddc==0; 
    formulastring='';
else
    formulastring=[num2str(fun.overaddc) ' + ']; 
end
if fun.overmultc~=1; 
    formulastring=[formulastring num2str(fun.overmultc) ' ']; 
end
if param.dim~=1; 
    formulastring=[formulastring 'prod_{j=1}^' int2str(param.dim)];
end
formulastring=[formulastring '['];
numdim=min(6,param.dim);
if param.dim<=numdim
    paramendstring=char(10); 
else
    paramendstring=[' ...' char(10)];
end
if fun.addc==0; 
    paramstring='';
else
    formulastring=[formulastring 'a_j + '];
    paramstring=['    a_j = ' num2str(fun.addc(1:numdim),'%6.3g')...
        paramendstring]; 
end
if fun.scale~=1; 
    formulastring=[formulastring 'b_j']; 
    paramstring=[paramstring ...
        '    b_j = ' num2str(fun.scale(1:numdim),'%6.3g')...
        paramendstring]; 
end
formulastring=[formulastring ...
    'Indicator_[0,p_j](x - z_j mod(up_j-lo_j)+lo_j'];
paramstring=[paramstring ...
    '    p_j = ' num2str(fun.shape(1:numdim),'%6.3g') paramendstring ...
    '    z_j = ' num2str(fun.shift(1:numdim),'%6.3g') paramendstring ...
    '    lo_j = ' num2str(param.interval(1,1:numdim),'%6.3g') ...
    paramendstring...
    '    hi_j = ' num2str(param.interval(2,1:numdim),'%6.3g') ...
    paramendstring]; 
param.funDescribe=['   f(x) = ' formulastring ']' char(10) paramstring];
end

function f=stepfun(x,overaddc,overmultc,addc,scale,shape,shift,bmina)
    n=size(x,1);
    f=overaddc+overmultc.*prod(repmat(addc,n,1) ...
        + repmat(scale,n,1)...
        .* (mod(x-repmat(shift,n,1),repmat(bmina,n,1))...
        <=repmat(shape,n,1)),2);
end
%% Gaussian Function
function [testfun,param]=makeGaussianTestFun(fun,param)
%Create the exponential test function
fun=verifyparam(fun,{'shape','scale','addc','overaddc','overmultc','center'}, ...
    {[1 param.dim],[1 param.dim],[1 param.dim],[1 1],[1 1],[1 param.dim]}, ...
    {1,1,0,0,1,0.5});
testfun=@(x) gaussianfun(x,fun.overaddc,fun.overmultc,fun.addc,fun.scale,...
    fun.shape,fun.center);

%% Compute exact integral of this function
if strcmp(param.measure,'uniform')
    bmina=param.interval(2,:)-param.interval(1,:);
    
    param.gaussian1=sqrt(pi)*fun.shape./(2*bmina).*...
        (erf((param.interval(2,:)-fun.center)...
        ./fun.shape)-erf((param.interval(1,:)-fun.center)./fun.shape));
    
    param.gaussian2=sqrt(pi)*fun.shape./(2*sqrt(2)*bmina).*...
        (erf(sqrt(2)*(param.interval(2,:)-fun.center)...
        ./fun.shape) - erf(sqrt(2)*(param.interval(1,:)-fun.center)./fun.shape));
    param.gaussian3=sqrt(pi)*fun.shape./(2*sqrt(3)*bmina).*...
        (erf(sqrt(3)*(param.interval(2,:)-fun.center)...
        ./fun.shape) - erf(sqrt(3)*(param.interval(1,:)-fun.center)./fun.shape));
    param.gaussian4=sqrt(pi)*fun.shape./(2*2*bmina).*...
        (erf(2*(param.interval(2,:)-fun.center)...
        ./fun.shape) - erf(2*(param.interval(1,:)-fun.center)./fun.shape));
            
    param.exactintegral=fun.overaddc.*prod(bmina)+ fun.overmultc.*prod(fun.addc ...
        .*bmina + fun.scale.*param.gaussian1.*bmina);
    
   param.exactmean=fun.overaddc+ fun.overmultc.*prod(fun.addc ...
        + fun.scale.*param.gaussian1);    
    
    param.moment1=fun.overmultc.*prod(fun.addc ...
        + fun.scale.*param.gaussian1);
    
    param.moment2=(fun.overmultc.^2.*...
        prod(fun.addc.^2+2*fun.addc.*fun.scale.*param.gaussian1...
        +fun.scale.^2.*param.gaussian2));
  
    param.exactvariance=param.moment2-param.moment1.^2;
    
    param.moment3=(fun.overmultc.^3.*prod(fun.addc.^3+...
        3*fun.addc.^2.*fun.scale.*param.gaussian1+...
        3*fun.addc.*fun.scale.^2.*param.gaussian2+...
        fun.scale.^3.*param.gaussian3));
    
    param.moment4=(fun.overmultc.^4.*prod(fun.addc.^4+...
        4*fun.addc.^3.*fun.scale.*param.gaussian1+...
        6*fun.addc.^2.*fun.scale.^2.*param.gaussian2+...
        4*fun.addc.*fun.scale.^3.*param.gaussian3+...
        fun.scale.^4.*param.gaussian4));
    
    param.exactkurtosis=(param.moment4-4*param.moment1...
        .*param.moment3+6*param.moment1.^2.*param.moment2...
        -4*param.moment1.^4+param.moment1.^4)...
        ./(param.exactvariance.^2);
    
    
else %normal distribution
    param.exactintegral=fun.overaddc ...
        + fun.overmultc.*prod(fun.addc ...
        + fun.scale.*exp(fun.shape.^2/2));
    exactintegralminus=fun.overmultc.*prod(fun.addc ...
        + fun.scale.*exp(fun.shape.^2/2));
    fminusint2=fun.overmultc.^2.*prod(fun.addc.^2 ...
        + 2.*fun.addc.*fun.scale.*exp(fun.shape.^2/2)...
        + fun.scale.^2.*exp((2*fun.shape).^2/2));
    fminusint3=fun.overmultc.^3.*prod(fun.addc.^3 ...
        + 3.*fun.addc.^2.*fun.scale.*exp(fun.shape.^2/2)...
        + 3*fun.addc.*fun.scale.^2.*exp((2*fun.shape).^2/2)...
        +fun.scale.^3.*exp((3*fun.shape).^2/2));
    fminusint4=fun.overmultc.^4.*prod(fun.addc.^4 ...
        + 4.*fun.addc.^3.*fun.scale.*exp(fun.shape.^2/2)...
        + 6*fun.addc.^2.*fun.scale.^2.*exp((2*fun.shape).^2/2)...
        +fun.addc.*fun.scale.^3.*exp((3*fun.shape).^2/2)...
        +fun.scale.^4.*exp(4*fun.shape).^2/2);
    param.exactvariance=fminusint2-exactintegralminus.^2;
      param.exactkurtosis=(fminusint4-4*exactintegralminus...
        .*fminusint3+6*exactintegralminus.^2.*fminusint2...
        -4*exactintegralminus.^4+exactintegralminus.^4)...
        ./(param.exactvariance.^2);
    
end

%% Prepare text description of the exponential function
if fun.overaddc==0; 
    formulastring='';
else
    formulastring=[num2str(fun.overaddc) ' + ']; 
end
if fun.overmultc~=1; 
    formulastring=[formulastring num2str(fun.overmultc) ' ']; 
end
if param.dim~=1; 
    formulastring=[formulastring 'prod_{j=1}^' int2str(param.dim)];
end
formulastring=[formulastring '['];
numdim=min(6,param.dim);
if param.dim<=numdim
    paramendstring=char(10); 
else
    paramendstring=[' ...' char(10)];
end
if fun.addc==0; 
    paramstring='';
else
    formulastring=[formulastring 'a_j + '];
    paramstring=['    a_j = ' num2str(fun.addc(1:numdim),'%6.3g')...
        paramendstring]; 
end
if fun.scale~=1; 
    formulastring=[formulastring 'b_j']; 
    paramstring=[paramstring ...
        '    b_j = ' num2str(fun.scale(1:numdim),'%6.3g')...
        paramendstring]; 
end
if fun.center~=0;
    formulastring=[formulastring 'H_j'];
    paramstring=[paramstring 'H_j = '...
        num2str(fun.center(1:numdim),'%6.3g ') paramendstring];
if fun.shape==1; 
    formulastring=[formulastring 'exp(-(x_j-H_j)^2)'];
else
    formulastring=[formulastring 'exp(-(x_j-H_j)^2/c_j^2)'];
    paramstring=[paramstring ...
        '    c_j = ' num2str(fun.shape(1:numdim),'%6.3g')...
        paramendstring]; 
end
param.funDescribe=['   f(x) = ' formulastring ']' char(10) paramstring];
end
end

function f=gaussianfun(x,overaddc,overmultc,addc,scale,shape,center)
    n=size(x,1);
    f=overaddc+overmultc.*prod(repmat(addc,n,1) ...
        + repmat(scale,n,1).*exp(-(x-repmat(center,n,1)).^2./repmat(shape.^2,n,1)),2);
end    


%% Gaussian Function without add parameter
function [testfun,param]=makeGaussianKerTestFun(fun,param)
%Create the exponential test function
fun=verifyparam(fun,{'shape','scale','addc','overmultc','center'}, ...
    {[1 param.dim],[1 param.dim],[1 param.dim],[1 1],[1 param.dim]}, ...
    {1,1,0,0,1,0.5});
testfun=@(x) gaussiankerfun(x,fun.overmultc,fun.addc,fun.scale,...
    fun.shape,fun.center);

%% Compute exact integral of this function
if strcmp(param.measure,'uniform')
    bmina=param.interval(2,:)-param.interval(1,:);
    
    param.gaussian1=sqrt(pi)*fun.shape./(2*bmina).*...
        (erf((param.interval(2,:)-fun.center)...
        ./fun.shape)-erf((param.interval(1,:)-fun.center)./fun.shape));
    
    param.gaussian2=sqrt(pi)*fun.shape./(2*sqrt(2)*bmina).*...
        (erf(sqrt(2)*(param.interval(2,:)-fun.center)...
        ./fun.shape) - erf(sqrt(2)*(param.interval(1,:)-fun.center)./fun.shape));
    param.gaussian3=sqrt(pi)*fun.shape./(2*sqrt(3)*bmina).*...
        (erf(sqrt(3)*(param.interval(2,:)-fun.center)...
        ./fun.shape) - erf(sqrt(3)*(param.interval(1,:)-fun.center)./fun.shape));
    param.gaussian4=sqrt(pi)*fun.shape./(2*2*bmina).*...
        (erf(2*(param.interval(2,:)-fun.center)...
        ./fun.shape) - erf(2*(param.interval(1,:)-fun.center)./fun.shape));
            
    param.exactintegral=fun.overmultc.*prod(fun.addc ...
        .*bmina + fun.scale.*param.gaussian1.*bmina);
    
   param.exactmean=fun.overmultc.*prod(fun.addc ...
        + fun.scale.*param.gaussian1);    
    
    param.moment1=fun.overmultc.*prod(fun.addc ...
        + fun.scale.*param.gaussian1);
    
    param.moment2=(fun.overmultc.^2.*...
        prod(fun.addc.^2+2*fun.addc.*fun.scale.*param.gaussian1...
        +fun.scale.^2.*param.gaussian2));
  
    param.exactvariance=param.moment2-param.moment1.^2;
    
    param.moment3=(fun.overmultc.^3.*prod(fun.addc.^3+...
        3*fun.addc.^2.*fun.scale.*param.gaussian1+...
        3*fun.addc.*fun.scale.^2.*param.gaussian2+...
        fun.scale.^3.*param.gaussian3));
    
    param.moment4=(fun.overmultc.^4.*prod(fun.addc.^4+...
        4*fun.addc.^3.*fun.scale.*param.gaussian1+...
        6*fun.addc.^2.*fun.scale.^2.*param.gaussian2+...
        4*fun.addc.*fun.scale.^3.*param.gaussian3+...
        fun.scale.^4.*param.gaussian4));
    
    param.exactkurtosis=(param.moment4-4*param.moment1...
        .*param.moment3+6*param.moment1.^2.*param.moment2...
        -4*param.moment1.^4+param.moment1.^4)...
        ./(param.exactvariance.^2);    
else %normal distribution
    param.exactintegral=fun.overmultc.*prod(fun.addc ...
        + fun.scale.*exp(fun.shape.^2/2));
    exactintegralminus=fun.overmultc.*prod(fun.addc ...
        + fun.scale.*exp(fun.shape.^2/2));
    fminusint2=fun.overmultc.^2.*prod(fun.addc.^2 ...
        + 2.*fun.addc.*fun.scale.*exp(fun.shape.^2/2)...
        + fun.scale.^2.*exp((2*fun.shape).^2/2));
    fminusint3=fun.overmultc.^3.*prod(fun.addc.^3 ...
        + 3.*fun.addc.^2.*fun.scale.*exp(fun.shape.^2/2)...
        + 3*fun.addc.*fun.scale.^2.*exp((2*fun.shape).^2/2)...
        +fun.scale.^3.*exp((3*fun.shape).^2/2));
    fminusint4=fun.overmultc.^4.*prod(fun.addc.^4 ...
        + 4.*fun.addc.^3.*fun.scale.*exp(fun.shape.^2/2)...
        + 6*fun.addc.^2.*fun.scale.^2.*exp((2*fun.shape).^2/2)...
        +fun.addc.*fun.scale.^3.*exp((3*fun.shape).^2/2)...
        +fun.scale.^4.*exp(4*fun.shape).^2/2);
    param.exactvariance=fminusint2-exactintegralminus.^2;
      param.exactkurtosis=(fminusint4-4*exactintegralminus...
        .*fminusint3+6*exactintegralminus.^2.*fminusint2...
        -4*exactintegralminus.^4+exactintegralminus.^4)...
        ./(param.exactvariance.^2);
    
end

%% Prepare text description of the Gaussian function
if fun.overmultc~=1; 
    formulastring=[num2str(fun.overmultc) ' ']; 
end
if param.dim~=1; 
    formulastring=[formulastring 'prod_{j=1}^' int2str(param.dim)];
end
formulastring=[formulastring '['];
numdim=min(6,param.dim);
if param.dim<=numdim
    paramendstring=char(10); 
else
    paramendstring=[' ...' char(10)];
end
if fun.addc==0; 
    paramstring='';
else
    formulastring=[formulastring 'a_j + '];
    paramstring=['    a_j = ' num2str(fun.addc(1:numdim),'%6.3g')...
        paramendstring]; 
end
if fun.scale~=1; 
    formulastring=[formulastring 'b_j']; 
    paramstring=[paramstring ...
        '    b_j = ' num2str(fun.scale(1:numdim),'%6.3g')...
        paramendstring]; 
end
if fun.center~=0;
    formulastring=[formulastring 'H_j'];
    paramstring=[paramstring 'H_j = '...
        num2str(fun.center(1:numdim),'%6.3g ') paramendstring];
if fun.shape==1; 
    formulastring=[formulastring 'exp(-(x_j-H_j)^2)'];
else
    formulastring=[formulastring 'exp(-(x_j-H_j)^2/c_j^2)'];
    paramstring=[paramstring ...
        '    c_j = ' num2str(fun.shape(1:numdim),'%6.3g')...
        paramendstring]; 
end
param.funDescribe=['   f(x) = ' formulastring ']' char(10) paramstring];
end
end

function f=gaussiankerfun(x,overmultc,addc,scale,shape,center)
    n=size(x,1);
    f=overmultc.*prod(repmat(addc,n,1) ...
        + repmat(scale,n,1).*exp(-(x-repmat(center,n,1)).^2./repmat(shape.^2,n,1)),2);
end    

    

