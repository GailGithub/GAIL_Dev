%Test the new cubMC routine
function [res,out_param] = TestcubMCDiffSettings(test,fun,in_param,out_param)
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
nsigold=in_param.n0;

for irep=1:test.nrep
    if round(irep/test.howoftenrep)==irep/test.howoftenrep, irep, end
    [testfun,fun,out_param]=test.randchoicefun(fun,in_param,test.randch,irep);
    if strcmp(fun.funtype,'step')
        res.exactkurtosis(irep)=out_param.exactkurtosis;
        res.exactvariance(irep)=out_param.exactvariance;
    end
    if strcmp(fun.funtype,'exp')
        res.exactkurtosis(irep)=out_param.exactkurtosis;
        res.exactvariance(irep)=out_param.exactvariance;

    end
    if strcmp(fun.funtype,'gaussian')
        res.exactkurtosis(irep)=out_param.exactkurtosis;
        res.exactvariance(irep)=out_param.exactvariance;

    end
    
    % Evaluate integral for iid
    if any(strcmp('iid',test.whichsample))
        in_param.n0=nsigold;
        %in_param.sample='iid';
        [~,out_param]=cubMC_g(testfun,in_param.interval,in_param);
        %keyboard;
        if irep==1; 
            res.iidkurtmax=out_param.kurtmax; 
        end
        res.iidexit(irep)=out_param.exit;
        res.iidQ(irep)=out_param.Q;
        res.iiderr(irep)=abs(out_param.exactintegral-out_param.Q);
        res.iidtime(irep)=out_param.time;
        res.iidneval(irep)=out_param.n;
        res.estvariance(irep)=out_param.estvari;

    end

    % Evaluate integral for heavy duty iid
    if any(strcmp('iidheavy',test.whichsample))
        in_param.n0=nsigold*2^5; %larger n to compute sigma
        %param.sample='iid';
        [~,out_param]=cubMC_g(testfun,in_param.interval,in_param);
        if irep==1; res.iidheavykurtmax=out_param.kurtmax; end
        res.iidheavyexit(irep)=out_param.exit;
        res.iidheavyQ(irep)=out_param.Q;
        res.iidheavyerr(irep)=abs(out_param.exactintegral-out_param.Q);
        res.iidheavytime(irep)=out_param.time;
        res.iidheavyneval(irep)=out_param.n;
        res.estvariance(irep)=out_param.estvari;

    end
    
%     % Evaluate integral for Sobol
%     if any(strcmp('Sobol',test.whichsample))
%         in_param.n0=nsigold;
%         in_param.sample='sobol';
%         in_param.scramble=true;
%         [~,out_param]=cubMC(testfun,in_param.interval,in_param);
%         res.Sobolexit(irep)=out_param.exit;
%         res.SobolQ(irep)=out_param.Q;
%         res.Sobolerr(irep)=abs(out_param.exactintegral-out_param.Q);
%         res.Soboltime(irep)=out_param.time;
%         res.Sobolneval(irep)=out_param.n;
% 
%     end
% 
%     % Evaluate integral for Sobol heavy
%     if any(strcmp('Sobolheavy',test.whichsample))
%         in_param.n0=nsigold*2^5; %larger n to compute sigma
%         in_param.sample='sobol';
%         in_param.scramble=true;
%         [~,out_param]=cubMC(testfun,in_param.interval,in_param);
%         res.Sobolheavyexit(irep)=out_param.exit;
%         res.SobolheavyQ(irep)=out_param.Q;
%         res.Sobolheavyerr(irep)=abs(out_param.exactintegral-out_param.Q);
%         res.Sobolheavytime(irep)=out_param.time;
%         res.Sobolheavyneval(irep)=out_param.n;
% 
%     end

    % Evaluate integral for MATLAB quad
    if in_param.dim==1
        if any(strcmp('quad',test.whichsample))
            tic,
            res.quadQ(irep)=quad(@(x) transpose(testfun(x')),...
                in_param.interval(1),in_param.interval(2),in_param.tol);
            res.quaderr(irep)=abs(out_param.exactintegral-res.quadQ(irep));
            res.quadtime(irep)=toc;
        end
        if any(strcmp('quadgk',test.whichsample))
            tic,
            res.quadgkQ(irep)=quadgk(@(x) transpose(testfun(x')),...
                in_param.interval(1),in_param.interval(2),'AbsTol',in_param.tol); 
            res.quadgkerr(irep)=abs(out_param.exactintegral-res.quadgkQ(irep));
            res.quadgktime(irep)=toc;
        end
    end
   % Evaluate integral for chebfun
    if in_param.dim==1
        if any(strcmp('chebfun',test.whichsample))
            tic,
            res.chebfunQ(irep)=sum(chebfun(@(x) testfun(x),...
                [in_param.interval(1),in_param.interval(2)],'minsamples',9,...
                'splitting','on'),in_param.interval(1),in_param.interval(2));
            res.chebfunerr(irep)=abs(out_param.exactintegral-res.chebfunQ(irep));
            res.chebfuntime(irep)=toc;
        end
        if any(strcmp('chebfunheavy',test.whichsample))
            tic,
            %chebfunpref('minsamples',1025)
            res.chebfunheavyQ(irep)=sum(chebfun(@(x) testfun(x),...
                [in_param.interval(1),in_param.interval(2)],'minsamples',257,...
                'splitting','on'),in_param.interval(1),in_param.interval(2));
            res.chebfunheavyerr(irep)=abs(out_param.exactintegral-res.chebfunheavyQ(irep));
            res.chebfunheavytime(irep)=toc;
        end        
    end    
end

timestamp=datestr(now);
timestamp(timestamp==' ')='_';
timestamp(timestamp==':')='.';
save(['./Results/TestcubMCon-' fun.funtype '-' in_param.measure '-Out-' ...
    timestamp 'N' int2str(test.nrep) 'd' int2str(in_param.dim)  ...
    'tol' num2str(in_param.tol) '.mat']) 

toc(tstartwhole)