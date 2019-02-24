%Test the new cubMC routine
function res=TestSettings_rel(test,fun,param)
tstartwhole=tic;

% Initialize variables
tempinitial=zeros(test.nrep,1);
res.dim=tempinitial;
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
if any(strcmp('cubSobol',test.whichsample))
    res.Sobolexit=tempinitial;
    res.SobolQ=tempinitial;
    res.Sobolerr=tempinitial;
    res.Soboltime=tempinitial;
    res.Sobolneval=tempinitial;
end
if any(strcmp('cubLattice',test.whichsample))
    res.Latticeexit=tempinitial;
    res.LatticeQ=tempinitial;
    res.Latticeerr=tempinitial;
    res.Latticetime=tempinitial;
    res.Latticeneval=tempinitial;
end

for irep=1:test.nrep
    if round(irep/test.howoftenrep)==irep/test.howoftenrep, irep, end
    %keyboard
    if any(strcmp(test.whichsample,{'iid','iidheavy'}))
       param.sample='iid';
       [testfuniid,fun,param]= ...
          test.randchoicefun(fun,param,test.randch,irep);
    elseif any(strcmp(test.whichsample,{'cubSobol','cubLattice'}))
       param.sample='qmc';
       [testfunqmc,fun,param]= ...
          test.randchoicefun(fun,param,test.randch,irep);
    else
       error('Don''t know the sample type')
    end
    res.dim(irep)=param.dim;
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
        param.n_sigma=1e4;
        param.sample='iid';
        [Q,out_param]=cubMC_g(testfuniid,param.interval,param);
        if irep==1; 
            res.iidkurtmax=out_param.kurtmax; 
        end
        res.iidexit(irep)=out_param.exit;
        res.iidQ(irep)=Q;
        res.iiderr(irep)=abs(param.exactintegral-Q);
        res.iidtime(irep)=out_param.time;
        res.iidneval(irep)=out_param.n;
    end

    % Evaluate integral for heavy duty iid
    if any(strcmp('iidheavy',test.whichsample))
        param.n_sigma=1e5; %larger n to compute sigma
        param.sample='iid';
        [Q,out_param]=cubMC_g(testfuniid,param.interval,param);
        if irep==1; res.iidheavykurtmax=out_param.kurtmax; end
        res.iidheavyexit(irep)=out_param.exit;
        res.iidheavyQ(irep)=Q;
        res.iidheavyerr(irep)=abs(param.exactintegral-Q);
        res.iidheavytime(irep)=out_param.time;
        res.iidheavyneval(irep)=out_param.n;
    end
    
    % Evaluate integral using cubSobol
    if any(strcmp('cubSobol',test.whichsample))
        [q,out_param]=...
           cubSobol_g(testfunqmc,[zeros(1,param.dim);ones(1,param.dim)],...
           'abstol',param.abstol,'reltol',param.reltol,'measure',param.measure,...
           'mmax',param.mmax);
%         res.Sobolexit(irep)=out_param.overbudget;
        res.SobolQ(irep)=q;
        res.Sobolexact = param.exactintegral;
        res.Sobolerr(irep)=abs(param.exactintegral-q);
        res.Soboltime(irep)=out_param.time;
        res.Sobolneval(irep)=out_param.n;
    end
    
    % Evaluate integral using cubSobol
    if any(strcmp('cubLattice',test.whichsample))
        [q,out_param]=...
           cubLattice_g(testfunqmc,[zeros(1,param.dim);ones(1,param.dim)],...
           'abstol',param.abstol,'reltol',param.reltol,'measure',param.measure,...
           'mmax',param.mmax,'transform',param.transform);
%         res.Latticeexit(irep)=out_param.overbudget;
        res.LatticeQ(irep)=q;
        res.Latticeexact = param.exactintegral;
        res.Latticeerr(irep)=abs(param.exactintegral-q);
        res.Latticetime(irep)=out_param.time;
        res.Latticeneval(irep)=out_param.n;
    end  
end

timestamp=datestr(now,'yyyy-mm-dd-HH-MM');
save(['TestCubature-' fun.funtype '-' param.measure ...
   '-N-' int2str(test.nrep)  ...
    '-' param.method '_rel.mat'])


toc(tstartwhole)