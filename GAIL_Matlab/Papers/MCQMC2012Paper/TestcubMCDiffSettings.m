%Test the new cubMC routine
function res=TestcubMCDiffSettings(test,fun,param)
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
if any(strcmp('iid',test.whichsample))
    res.iidexit=tempinitial;
    res.iidQ=tempinitial;
    res.iiderr=tempinitial;
    res.iidtime=tempinitial;
    res.iidneval=tempinitial;
end
if any(strcmp('iidheavy',test.whichsample))
    res.iidheavyexit=tempinitial;
    res.iidheavyQ=tempinitial;
    res.iidheavyerr=tempinitial;
    res.iidheavytime=tempinitial;
    res.iidheavyneval=tempinitial;
end
if any(strcmp('Sobol',test.whichsample))
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
if any(strcmp('quad',test.whichsample))
    res.quadQ=tempinitial;
    res.quaderr=tempinitial;
    res.quadtime=tempinitial;
end
if any(strcmp('quadgk',test.whichsample))
    res.quadgkQ=tempinitial;
    res.quadgkerr=tempinitial;
    res.quadgktime=tempinitial;
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
nsigold=param.n0;

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
    
    % Evaluate integral for iid
    if any(strcmp('iid',test.whichsample))
        param.n0=nsigold;
        param.sample='iid';
        [~,param]=cubMC(testfun,param.interval,param);
        if irep==1; 
            res.iidkurtmax=param.kurtmax; 
        end
        res.iidexit(irep)=param.exit;
        res.iidQ(irep)=param.Q;
        res.iiderr(irep)=abs(param.exactintegral-param.Q);
        res.iidtime(irep)=param.time;
        res.iidneval(irep)=param.n;
        res.estvariance(irep)=param.estvari;

    end

    % Evaluate integral for heavy duty iid
    if any(strcmp('iidheavy',test.whichsample))
        param.n0=nsigold*2^5; %larger n to compute sigma
        param.sample='iid';
        [~,param]=cubMC(testfun,param.interval,param);
        if irep==1; res.iidheavykurtmax=param.kurtmax; end
        res.iidheavyexit(irep)=param.exit;
        res.iidheavyQ(irep)=param.Q;
        res.iidheavyerr(irep)=abs(param.exactintegral-param.Q);
        res.iidheavytime(irep)=param.time;
        res.iidheavyneval(irep)=param.n;
        res.estvariance(irep)=param.estvari;

    end
    
    % Evaluate integral for Sobol
    if any(strcmp('Sobol',test.whichsample))
        param.n0=nsigold;
        param.sample='sobol';
        param.scramble=true;
        [~,param]=cubMC(testfun,param.interval,param);
        res.Sobolexit(irep)=param.exit;
        res.SobolQ(irep)=param.Q;
        res.Sobolerr(irep)=abs(param.exactintegral-param.Q);
        res.Soboltime(irep)=param.time;
        res.Sobolneval(irep)=param.n;

    end

    % Evaluate integral for Sobol heavy
    if any(strcmp('Sobolheavy',test.whichsample))
        param.n0=nsigold*2^5; %larger n to compute sigma
        param.sample='sobol';
        param.scramble=true;
        [~,param]=cubMC(testfun,param.interval,param);
        res.Sobolheavyexit(irep)=param.exit;
        res.SobolheavyQ(irep)=param.Q;
        res.Sobolheavyerr(irep)=abs(param.exactintegral-param.Q);
        res.Sobolheavytime(irep)=param.time;
        res.Sobolheavyneval(irep)=param.n;

    end

    % Evaluate integral for MATLAB quad
    if param.dim==1
        if any(strcmp('quad',test.whichsample))
            tic,
            res.quadQ(irep)=quad(@(x) transpose(testfun(x')),...
                param.interval(1),param.interval(2),param.tol);
            res.quaderr(irep)=abs(param.exactintegral-res.quadQ(irep));
            res.quadtime(irep)=toc;
        end
        if any(strcmp('quadgk',test.whichsample))
            tic,
            res.quadgkQ(irep)=quadgk(@(x) transpose(testfun(x')),...
                param.interval(1),param.interval(2),'AbsTol',param.tol); 
            res.quadgkerr(irep)=abs(param.exactintegral-res.quadgkQ(irep));
            res.quadgktime(irep)=toc;
        end
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
subdir = 'MCQMC2012PaperOutput';
filename = ['TestcubMCon-' fun.funtype '-' param.measure  '-N'...
    int2str(test.nrep) 'd' int2str(param.dim)  ...
    'tol' num2str(param.tol)];
gail.save_mat(subdir, filename, true, res,test,fun,param);
toc(tstartwhole)