classdef meanMC
    
    properties (SetAccess=public)
        in_param = struct( ...
            'abstol', 1e-2, ... % default absolute error tolerance
            'reltol', 1e-1, ... % default relative error tolerance
            'alpha', 0.01, ... % default uncertainty
            'fudge', 1.2, ... % default standard deviation inflation factor
            'nb', 1e3, ... % initial sample size (per control variate)
            ... %for estimating the coefficients of control variates
            'nSig', 1e4, ... % default the sample size to estimate the variance
            'n1', 1e4, ... % default the initial sample size to estimate the mean
            'tbudget', 100, ... % default time budget
            'nbudget', 1e9) % default sample budget;
        method = 'na'
        Yrand = @(n) rand(n,1).^2
        YXrand = @(n) rand(n,1).^2
        muX = []
    end
    
    properties (Constant, Hidden) %do not change & not seen
        allowVarRed = {'no', 'cv'}
    end
    
    methods
        % Creating an meanMC process
        function obj = meanMC(varargin)
            if nargin > 0
                val = varargin{1};
                if isa(val,'meanMC')
                    obj.in_param = val.in_param;
                    obj.method = val.method;
                    obj.Yrand = val.Yrand;
                    obj.YXrand = val.YXrand;
                    obj.muX = val.muX;
                    if nargin == 1
                        return
                    else
                        val = varargin{2};
                    end
                end
                if isstruct(val)
                    if isfield(val,'in_param')
                        obj.in_param = val.in_param;
                    end
                    if isfield(val,'method')
                        obj.method = val.method;
                    end
                    if isfield(val,'Yrand')
                        obj.Yrand = val.Yrand;
                    end
                    if isfield(val,'YXrand')
                        obj.YXrand = val.YXrand;
                    end
                    if isfield(val,'muX')
                        obj.muX = val.muX;
                    end
                end
            end
            
        end
        
        % Set the properties of the meanMC object
        function obj = set.in_param(obj,val)
            if isfield(val,'abstol')
                validateattributes(val.abstol,{'numeric'}, ...
                    {'nonnegative'})
                obj.in_param.abstol =val.abstol; %row
            end
            if isfield(val,'reltol')
                validateattributes(val.reltol,{'numeric'}, ...
                    {'nonnegative'})
                assert(val.reltol <=1)
                obj.in_param.reltol =val.reltol; %row
            end
            if isfield(val,'alpha')
                validateattributes(val.alpha,{'numeric'}, ...
                    {'nonnegative'})
                assert(val.alpha <1)
                obj.in_param.alpha =val.alpha; %row
            end
        end
        
        function [tmu, out_param] = genMu(obj)
            [tmu, out_param]=meanMC_g(obj.Yrand, obj.in_param);
        end 
        
    end
    
    
    methods (Access = protected)
        function [tmu,out_param]=meanMC_g(Yrand, in_param)
            tstart = tic;
            out_param = in_param;
            n1 = 2;
            Yrand(n1); %let it run once to load all the data. warm up the machine.
            nsofar = n1;
            
            ntry = 10;% try several samples to get the time
            tic;
            Yrand(ntry);
            ttry=toc;
            tpern = ttry/ntry; % calculate time per sample
            nsofar = nsofar+ntry; % update n so far
            out_param.exit = 0;
            if tpern<1e-7;%each sample use very very little time
                booster = 8;
                tic;Yrand(ntry*booster);ttry2 = toc;
                ntry = ntry*[1 booster];
                ttry = [ttry ttry2];% take eight times more samples to try
            elseif tpern>=1e-7 && tpern<1e-5 %each sample use very little time
                booster = 6;
                tic;Yrand(ntry*booster);ttry2 = toc;
                ntry = ntry*[1 booster];
                ttry = [ttry ttry2];% take six times more samples to try
            elseif tpern>=1e-5 && tpern<1e-3 %each sample use little time
                booster = 4;
                tic;Yrand(ntry*booster);ttry2 = toc;
                ntry = ntry*[1 booster];
                ttry = [ttry ttry2];% take four times more samples to try
            elseif  tpern>=1e-3 && tpern<1e-1 %each sample use moderate time
                booster = 2;
                tic;Yrand(ntry*booster);ttry2 = toc;
                ntry = ntry*[1 booster];
                ttry = [ttry ttry2];% take two times more samples to try
            else %each sample use lots of time, stop try
            end
            [tmu,out_param] = meanmctolfun(Yrand,out_param,ntry,ttry,nsofar,tstart);
        end
        
        function [tmu,out_param] =  meanmctolfun(Yrand,out_param,ntry,ttry,nsofar,tstart)
            tic;
            Yval = Yrand(out_param.nSig);% get samples to estimate variance
            t_sig = toc;%get the time for calculating nSig function values.
            nsofar = nsofar+out_param.nSig;% update the samples that have been used
            out_param.nremain = gail.estsamplebudget(out_param.tbudget,...
                out_param.nbudget,[ntry out_param.nSig],nsofar,tstart,[ttry t_sig]);
            %update the nremain could afford until now
            out_param.var = var(Yval);% calculate the sample variance--stage 1
            sig0 = sqrt(out_param.var);% standard deviation
            sig0up = out_param.fudge*sig0;% upper bound on the standard deviation
            alpha_sig = out_param.alpha/2;% the uncertainty for variance estimation
            out_param.kurtmax = (out_param.nSig-3)/(out_param.nSig-1) ...
                + ((alpha_sig*out_param.nSig)/(1-alpha_sig))...
                *(1-1/out_param.fudge^2)^2;
            %the upper bound on the modified kurtosis
            npcmax = 1e6;%constant to do iteration and mean calculation4
            if out_param.reltol ==0
                out_param.tau = 1;
                alphai = 1-(1-out_param.alpha)/(1-alpha_sig);
                if sig0up == 0; % if the variance is zero, just take n_sigma samples
                    out_param.n = out_param.nSig;
                else
                    toloversig = out_param.abstol/sig0up;
                    % absolute error tolerance over sigma
                    out_param.n = nchebe(toloversig,alphai,out_param.kurtmax);
                    if out_param.n > out_param.nremain;
                        out_param.exit=1; %pass a flag
                        meanMC_g_err(out_param); % print warning message
                        out_param.n = out_param.nremain;% update n
                    end
                end
                tmu = gail.evalmean(Yrand,out_param.n,npcmax);%evaluate the mean
                nsofar = nsofar+out_param.n;%total sample used
            else
                alphai = (out_param.alpha-alpha_sig)/(2*(1-alpha_sig));
                %uncertainty to do iteration
                eps1 = ncbinv(out_param.n1,alphai,out_param.kurtmax);
                %tolerance for initial estimation
                out_param.tol(1) = sig0up*eps1;
                %the width of initial confidence interval for the mean
                i=1;
                out_param.n(i) = out_param.n1;% initial sample size to do iteration
                while true
                    out_param.tau = i;%step of the iteration
                    if out_param.n(i) > out_param.nremain;
                        % if the sample size used for initial estimation is
                        % larger than nremain, print warning message and use nremain
                        out_param.exit=1; %pass a flag
                        meanMC_g_err(out_param); % print warning message
                        out_param.n(i) = out_param.nremain;% update n
                        tmu = gail.evalmean(Yrand,out_param.n(i),npcmax);%evaluate the mean
                        nsofar = nsofar+out_param.n(i);%total sample used
                        break;
                    end
                    out_param.hmu(i) = gail.evalmean(Yrand,out_param.n(i),npcmax);%evaluate mean
                    nsofar = nsofar+out_param.n(i);
                    out_param.nremain = out_param.nremain-out_param.n(i);%update n so far and nremain
                    errtype = 'max';
                    % error type, see the function 'tolfun' at Algoithms/+gail/ directory
                    % for more info
                    theta  = 0;% relative error case
                    deltaplus = (gail.tolfun(out_param.abstol,out_param.reltol,...
                        theta,out_param.hmu(i) - out_param.tol(i),errtype)...
                        +gail.tolfun(out_param.abstol,out_param.reltol,...
                        theta,out_param.hmu(i) + out_param.tol(i),errtype))/2;
                    % a combination of tolfun, which used to decide stopping time
                    if deltaplus >= out_param.tol(i) % stopping criterion
                        deltaminus= (gail.tolfun(out_param.abstol,out_param.reltol,...
                            theta,out_param.hmu(i) - out_param.tol(i),errtype)...
                            -gail.tolfun(out_param.abstol,out_param.reltol,...
                            theta,out_param.hmu(i) + out_param.tol(i),errtype))/2;
                        % the other combination of tolfun, which adjust the hmu a bit
                        tmu = out_param.hmu(i)+deltaminus;
                        break;
                    else
                        out_param.tol(i+1) = min(out_param.tol(i)/2,max(out_param.abstol,...
                            0.95*out_param.reltol*abs(out_param.hmu(i))));
                        i=i+1;
                    end
                    toloversig = out_param.tol(i)/sig0up;%next tolerance over sigma
                    alphai = (out_param.alpha-alpha_sig)/(1-alpha_sig)*2.^(-i);
                    %update the next uncertainty
                    out_param.n(i) = nchebe(toloversig,alphai,out_param.kurtmax);
                end
            end
            %get the next sample size needed
            out_param.ntot = nsofar;%total sample size used
            out_param.time=toc(tstart); %elapsed time
        end
        
        function ncb = nchebe(toloversig,alpha,kurtmax)
            %this function uses Chebyshev and Berry-Esseen Inequality to calculate the
            %sample size needed
            ncheb = ceil(1/(toloversig^2*alpha));%sample size by Chebyshev's Inequality
            A=18.1139;
            A1=0.3328;
            A2=0.429; % three constants in Berry-Esseen inequality
            M3upper = kurtmax^(3/4);
            %the upper bound on the third moment by Jensen's inequality
            BEfun2=@(logsqrtn)gail.stdnormcdf(-exp(logsqrtn).*toloversig)...
                +exp(-logsqrtn).*min(A1*(M3upper+A2), ...
                A*M3upper./(1+(exp(logsqrtn).*toloversig).^3))-alpha/2;
            % Berry-Esseen function, whose solution is the sample size needed
            logsqrtnCLT=log(gail.stdnorminv(1-alpha/2)/toloversig);%sample size by CLT
            nbe=ceil(exp(2*fzero(BEfun2,logsqrtnCLT)));
            %calculate Berry-Esseen n by fzero function
            ncb = min(ncheb,nbe);%take the min of two sample sizes.
        end
        
        function eps = ncbinv(n1,alpha1,kurtmax)
            %This function calculate the reliable upper bound on error when given
            %Chebyshev and Berry-Esseen inequality and sample size n.
            NCheb_inv = 1/sqrt(n1*alpha1);
            % use Chebyshev inequality
            A=18.1139;
            A1=0.3328;
            A2=0.429; % three constants in Berry-Esseen inequality
            M3upper=kurtmax^(3/4);
            %using Jensen's inequality to bound the third moment
            BEfun=@(logsqrtb)gail.stdnormcdf(n1.*logsqrtb)...
                +min(A1*(M3upper+A2), ...
                A*M3upper./(1+(sqrt(n1).*logsqrtb).^3))/sqrt(n1)...
                - alpha1/2;
            % Berry-Esseen inequality
            logsqrtb_clt=log(sqrt(gail.stdnorminv(1-alpha1/2)/sqrt(n1)));
            %use CLT to get tolerance
            NBE_inv = exp(2*fzero(BEfun,logsqrtb_clt));
            %use fzero to get Berry-Esseen tolerance
            eps = min(NCheb_inv,NBE_inv);
            %take the min of Chebyshev and Berry Esseen tolerance
        end
        
    end
end
