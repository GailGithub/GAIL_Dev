%Test the new cubMC routine
function res=TestcubMCgDiffSettings(test,fun,param)
tstartwhole=tic;

% Initialize variables
tempinitial=zeros(test.nrep,1);
if strcmp(fun.funtype,'step')
    res.exactkurtosis=tempinitial;
    res.estvariance=tempinitial;
    res.exactvariance=tempinitial;
end
if strcmp(fun.funtype,'exp')
    res.exactkurtosis=tempinitial;
    res.estvariance=tempinitial;
    res.exactvariance=tempinitial;
    
end
if strcmp(fun.funtype,'gaussian')
    res.exactkurtosis=tempinitial;
    res.estvariance=tempinitial;
    res.exactvariance=tempinitial;
end

if strcmp(fun.funtype,'Keister')
    res.exactintegral = tempinitial;
end

if strcmp(fun.funtype,'product')
    res.exactkurtosis=tempinitial;
    res.estvariance=tempinitial;
    res.exactvariance=tempinitial;
end
if any(strcmp('iid',test.whichsample))
    res.iidexit=tempinitial;
    res.iidQ=tempinitial;
    res.iiderr=tempinitial;
    res.iidtime=tempinitial;
    res.iidneval=tempinitial;
    res.iidtau=tempinitial;
end
if any(strcmp('iidheavy',test.whichsample))
    res.iidheavyexit=tempinitial;
    res.iidheavyQ=tempinitial;
    res.iidheavyerr=tempinitial;
    res.iidheavytime=tempinitial;
    res.iidheavyneval=tempinitial;
    res.iidheavytau=tempinitial;
end
if any(strcmp('cubSobol',test.whichsample))
    res.Sobolexit=tempinitial;
    res.SobolQ=tempinitial;
    res.Sobolerr=tempinitial;
    res.Soboltime=tempinitial;
    res.Sobolneval=tempinitial;
end
if any(strcmp('Sobolheavy',test.whichsample))
    res.Sobolheavyexit=tempinitial;
    res.SobolheavyQ=tempinitial;
    res.Sobolheavyerr=tempinitial;
    res.Sobolheavytime=tempinitial;
    res.Sobolheavyneval=tempinitial;
end
if any(strcmp('cubLattice',test.whichsample))
    res.Latticeexit=tempinitial;
    res.LatticeQ=tempinitial;
    res.Latticeerr=tempinitial;
    res.Latticetime=tempinitial;
    res.Latticeneval=tempinitial;
end

if any(strcmp('integral',test.whichsample))
    res.integralQ=tempinitial;
    res.integralerr=tempinitial;
    res.integraltime=tempinitial;
end
if any(strcmp('integral_g',test.whichsample))
    res.integralgQ=tempinitial;
    res.integralgerr=tempinitial;
    res.integralgtime=tempinitial;
end
if any(strcmp('chebfun',test.whichsample))
    res.chebfunQ=tempinitial;
    res.chebfunerr=tempinitial;
    res.chebfuntime=tempinitial;
end
if any(strcmp('chebfunheavy',test.whichsample))
    res.chebfunheavyQ=tempinitial;
    res.chebfunheavyerr=tempinitial;
    res.chebfunheavytime=tempinitial;
end

for irep=1:test.nrep
    gail.print_iterations(irep, 'irep', true);
    [testfun,fun,param]=test.randchoicefun(fun,param,test.randch,irep);
    if strcmp(fun.funtype,'step')
        res.exactkurtosis(irep)=param.exactkurtosis;
        res.exactvariance(irep)=param.exactvariance;
    end
    if strcmp(fun.funtype,'exp')
        res.exactkurtosis(irep)=param.exactkurtosis;
        res.exactvariance(irep)=param.exactvariance;
        
    end
    if strcmp(fun.funtype,'gaussian')
        res.exactkurtosis(irep)=param.exactkurtosis;
        res.exactvariance(irep)=param.exactvariance;
    end
    
    if strcmp(fun.funtype,'Keister')
        res.exactintegral(irep) = param.exactintegral;
    end
    
    if strcmp(fun.funtype,'product')
        res.exactkurtosis(irep)=param.exactkurtosis;
        res.exactvariance(irep)=param.exactvariance;
        
    end
    
    if any(strcmp('iid',test.whichsample))
        param.nSig=2^13;
        param.sample='iid';
        [Q,out_param]=cubMC_g(testfun,param.interval,param);
        if irep==1;
            res.iidkurtmax=out_param.kurtmax;
        end
        res.iidexit(irep)=out_param.exit;
        res.iidQ(irep)=Q;
        res.iiderr(irep)=abs(param.exactintegral-Q);
        res.iidtime(irep)=out_param.time;
        res.iidneval(irep)=out_param.ntot;
        res.iidtau(irep) = out_param.tau;
    end
    
    % Evaluate integral for heavy duty iid
    if any(strcmp('iidheavy',test.whichsample))
        param.nSig=2^18; %larger n to compute sigma
        param.sample='iid';
        [Q,out_param]=cubMC_g(testfun,param.interval,param);
        if irep==1; res.iidheavykurtmax=out_param.kurtmax; end
        res.iidheavyexit(irep)=out_param.exit;
        res.iidheavyQ(irep)=Q;
        res.iidheavyerr(irep)=abs(param.exactintegral-Q);
        res.iidheavytime(irep)=out_param.time;
        res.iidheavyneval(irep)=out_param.ntot;
        res.iidheavytau(irep) = out_param.tau;
        
    end
    
    % Evaluate integral using cubSobol
    if any(strcmp('cubSobol',test.whichsample))
        [q,out_param]=...
            cubSobol_g(testfun,param.interval,param);
        res.Sobolexit(irep)=out_param.exitflag(1);
        res.SobolQ(irep)=q;
        res.Sobolerr(irep)=abs(param.exactintegral-q);
        res.Soboltime(irep)=out_param.time;
        res.Sobolneval(irep)=out_param.n;
    end
    %
    %Evaluate integral using cubSobol
    if any(strcmp('cubLattice',test.whichsample))
        [q,out_param]=...
            cubLattice_g(testfun,param.interval,param);
        %            cubLattice_g(testfunqmc,param.dim,...
        %            'abstol',param.abstol,'density',param.measure,'transform','Baker',...
        %            'fudge',@(x) 20*2^-x);
        res.Latticeexit(irep)=out_param.exitflag(1);
        res.LatticeQ(irep)=q;
        res.Latticeerr(irep)=abs(param.exactintegral-q);
        res.Latticetime(irep)=out_param.time;
        res.Latticeneval(irep)=out_param.n;
    end
    %
    % Evaluate integral for MATLAB quad
    if param.dim==1
        if any(strcmp('integral',test.whichsample))
            tic,
            res.integralQ(irep)=integral(@(x) transpose(testfun(x')),...
                param.interval(1),param.interval(2),'AbsTol',param.abstol);
            res.integralerr(irep)=abs(param.exactintegral-res.integralQ(irep));
            res.integraltime(irep)=toc;
        end
        %         if any(strcmp('quadgk',test.whichsample))
        %             tic,
        %             res.quadgkQ(irep)=quadgk(@(x) transpose(testfun(x')),...
        %                 param.interval(1),param.interval(2),'AbsTol',param.abstol);
        %             res.quadgkerr(irep)=abs(param.exactintegral-res.quadgkQ(irep));
        %             res.quadgktime(irep)=toc;
        %         end
        %         if any(strcmp('integral_g',test.whichsample))
        %             tic,
        %             res.integralgQ(irep)=integral_g(@(x)transpose(testfun(x)),...
        %                 param.interval(1),param.interval(2),param.abstol);
        %             res.integralgerr(irep)=abs(param.exactintegral-res.integralgQ(irep));
        %             res.integralgtime(irep)=toc;
        %         end
    end
    % Evaluate integral for chebfun
    if param.dim==1
        if any(strcmp('chebfun',test.whichsample))
            tic,
            res.chebfunQ(irep)=sum(chebfun(@(x) testfun(x),...
                [param.interval(1),param.interval(2)],'minsamples',9,...
                'splitting','on'),param.interval(1),param.interval(2));
            res.chebfunerr(irep)=abs(param.exactintegral-res.chebfunQ(irep));
            res.chebfuntime(irep)=toc;
        end
        if any(strcmp('chebfunheavy',test.whichsample))
            tic,
            %chebfunpref('minsamples',1025)
            res.chebfunheavyQ(irep)=sum(chebfun(@(x) testfun(x),...
                [param.interval(1),param.interval(2)],'minsamples',257,...
                'splitting','on'),param.interval(1),param.interval(2));
            res.chebfunheavyerr(irep)=abs(param.exactintegral-res.chebfunheavyQ(irep));
            res.chebfunheavytime(irep)=toc;
        end
    end
end
disp(' ')
subdir = 'LanThesisOutput';
filename = ['TestcubMCon-' fun.funtype '-' param.measure  '-N'...
    int2str(test.nrep) 'd' int2str(param.dim)  ...
    'abstol' num2str(param.abstol)];
gail.save_mat(subdir, filename, true, res,test,fun,param);
toc(tstartwhole)