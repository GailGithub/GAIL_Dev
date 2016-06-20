classdef meanMC
    
    %% meanMC
    % is a class that uses Monte Carlo method to estimate the mean of a
    % random variable.
    %
    % Example 1
    % >> obj = meanMC
    % obj =
    %      meanMC with properties:
    %
    %             in_param_abstol: 1e-2
    %             in_param_abstol: 1e-1
    %              in_param_alpha: 0.01
    %              in_param_fudge: 1.2
    %               in_param_nSig: 1e4
    %                 in_param_n1: 1e4
    %            in_param_tbudget: 100
    %            in_param_nbudget: 1e9
    %                      method: {'plain'}
    %                          nc: 2e3
    %                       Yrand: @(n)rand(n,1).^2
    %             cv_param_YXrand: @(n)rand(n,1).^2
    %                cv_param_muX: []
    %                cv_param_ncv: 1e3
    %              cv_param_ridge: 0
    
    % Authors: Tianpei Qian, ... (and authors for the function |meanMC_g|)
    
    %% Properties
    % Below are the proporties for this class, along with their default
    % values
    
    properties (SetAccess=public)
        in_param = struct( ...
            'abstol', 1e-2, ... % default absolute error tolerance
            'reltol', 1e-1, ... % default relative error tolerance
            'alpha', 0.01, ... % default uncertainty
            'fudge', 1.2, ... % default standard deviation inflation factor
            'nSig', 1e4, ... % default sample size to estimate the variance
            'n1', 1e4, ... % default initial sample size to estimate the mean
            'tbudget', 100, ... % default time budget
            'nbudget', 1e9) % default sample budget;  
        
        method = {'plain'}
        
        nc = 2e3 % number of samples used to compare different methods
        % (per method)
        
        Yrand = @(n) rand(n,1).^2 % function that generate samples of the
        % random variable  
        
        cv_param = struct(...
            'YXrand', @(n) rand(n,1).^2, ... % function that generate samples of the
            ... % random variable and control variates
            'muX', [], ... % mean of control variates
            'ncv', 1e3, ... % sample size (per control variate)
            ... %for estimating the coefficients of control variates
            'ridge', 0) % parameter of ridge regression
        
        av_param = struct(...
            'YYrand', @(n) rand(n,2).^2) % function that generate antithetic pairs
        
    end
    
    properties (Constant, Hidden) %do not change & not seen
        allowMethod = {'plain', 'cv', 'av'}
    end
    
    %% Methods
    % Main feature: the function |genMu| estimates the mean of the
    % user-specified variable.
    
    methods
        % Creating an meanMC process
        function obj = meanMC(varargin)
            if nargin > 0
                val = varargin{1};
                if isa(val,'meanMC')
                    obj.in_param = val.in_param;
                    obj.method = val.method;
                    obj.Yrand = val.Yrand;
                    obj.cv_param = val.cv_param;
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
                    if isfield(val,'nc')
                        obj.nc = val.nc;
                    end
                    if isfield(val,'Yrand')
                        obj.Yrand = val.Yrand;
                    end
                    if isfield(val,'cv_param')
                        obj.cv_param = val.cv_param;
                    end
                    if isfield(val,'av_param')
                        obj.av_param = val.av_param;
                    end
                end
            end
            
        end
        
        % Set in_param of the meanMC object
        function obj = set.in_param(obj,val)
            if isfield(val,'abstol')
                validateattributes(val.abstol,{'numeric'}, ...
                    {'nonnegative'}) % abstol > 0
                obj.in_param.abstol =val.abstol;
            end
            if isfield(val,'reltol')
                validateattributes(val.reltol,{'numeric'}, ...
                    {'nonnegative', '<=', 1}) % 0 <= reltol <=1
                obj.in_param.reltol =val.reltol;
            end
            if isfield(val,'alpha')
                validateattributes(val.alpha,{'numeric'}, ...
                    {'nonnegative','<',1}) % 0<= alpha < 1
                obj.in_param.alpha =val.alpha;
            end
            if isfield(val,'fudge')
                validateattributes(val.fudge,{'numeric'}, ...
                    {'nonnegative','>',1}) % fudge > 1
                obj.in_param.fudge =val.fudge;
            end
            if isfield(val,'nSig')
                validateattributes(val.nSig,{'numeric'}, ...
                    {'positive'})
                assert(gail.isposge30(val.nSig)) % nSig >= 30
                obj.in_param.nSig =val.nSig;
            end
            if isfield(val,'n1')
                validateattributes(val.n1,{'numeric'}, ...
                    {'positive'})
                assert(gail.isposge30(val.n1)) % n1 >= 30
                obj.in_param.n1 =val.n1;
            end
            if isfield(val,'tbudget')
                validateattributes(val.tbudget,{'numeric'}, ...
                    {'positive'}) % tbudget > 0
                obj.in_param.tbudget =val.tbudget;
            end
            if isfield(val,'nbudget')
                validateattributes(val.nbudget,{'numeric'}, ...
                    {'positive'})
                assert(gail.isposge30(val.nbudget)) % mbudget >= 30
                obj.in_param.nbudget =val.nbudget;
            end
        end
        
        
        % Set nc of the meanMC object
        function obj = set.nc(obj,val)
            validateattributes(val,{'numeric'}, ...
                {'positive'}) % nc >= 30
            assert(gail.isposge30(val))
            obj.nc=val;
        end
        
        % Set method of the meanMC object
        function obj = set.method(obj,val)
            assert(all(any(strcmp(repmat(val,numel(obj.allowMethod),1), ...
                repmat(obj.allowMethod',1,numel(val))),1),2))
            obj.method=val;
        end
        
        % Set Yrand of the meanMC object
        function obj = set.Yrand(obj,val)
            assert(gail.isfcn(val))
            validateattributes(val(5), {'numeric'}, ...
                {'size', [5,1]})
            obj.Yrand=val;
        end
        
        % Set cv_param of the meanMC object
        function obj = set.cv_param(obj,val)
            if isfield(val,'ncv')
                validateattributes(val.ncv,{'numeric'}, ...
                    {'nonnegative'})
                assert(gail.isposge30(val.ncv)) % ncv >= 30
                obj.cv_param.ncv = val.ncv;
            end
            if isfield(val, 'YXrand')
                assert(gail.isfcn(val.YXrand))
                validateattributes(val.YXrand(5), {'numeric'}, ...
                    {'nrows',5})
                obj.cv_param.YXrand = val.YXrand;
            end
            if isfield(val, 'muX')
                assert(isnumeric(val.muX))
                if numel(val.muX)>0
                validateattributes(val.muX, {'numeric'}, ...
                    {'nrows',1})
                end
                obj.cv_param.muX = val.muX;
            end
            assert(size(obj.cv_param.YXrand(5),2)-1 ...
                ==length(obj.cv_param.muX))
            % number of control variates must equal length of muX
            if isfield(val, 'ridge')
                assert(isnumeric(val.ridge))
                validateattributes(val.ridge, {'numeric'}, ...
                    {'nonnegative'})
                obj.cv_param.ridge = val.ridge;
            end
        end
        
        % Set av_param of the meanMC object
        function obj = set.av_param(obj,val)
            assert(gail.isfcn(val.YYrand))
            validateattributes(val.YYrand(5), {'numeric'}, ...
                {'size', [5,2]})
            obj.av_param.YYrand=val.YYrand;
        end
        
        % estimate mu
        function [tmu, out_param] = genMu(obj)
            tstart = tic; % start the clock
            
            plain = any(strcmp(obj.method, 'plain'));
            cv = any(strcmp(obj.method, 'cv')) ...
                && size(obj.cv_param.YXrand(1),2) > 1; % at least one control variate
            av = any(strcmp(obj.method, 'av'));
            
            nmethods = plain+cv*length(obj.cv_param.ridge)+av; % number of methods
            Yrand_all = cell(1, nmethods); % to contain functions asssociated with each method
            comparison = (nmethods > 1); % whether they are multiple methods
            
            index = 1; % starting index of Yrand_all
            opt = 1; % index of the optiomal method
            nextra = 0; % to count samples used but not counted by meanMC_g   
            
            if plain
                Yrand_all{index} = obj.Yrand;
                index = index+1;
            end
            
            if cv
                ridge = obj.cv_param.ridge;
                for i = 1:length(ridge)
                    Yrand_all{index} = genYrand_cv(obj, ridge(i));
                    index = index + 1;
                end
                nextra = nextra + size(obj.cv_param.YXrand(1),2) ...
                    * obj.cv_param.ncv * length(ridge);
            end
            
            if av
                Yrand_all{index} = genYrand_av(obj);
            end 
            
            if comparison
                [Yrand_opt, index_opt] = selectYrand(obj,Yrand_all); % select the best method
                nextra = nextra + obj.nc * nmethods;
                index_list = 1:index;
                opt = index_list(index_opt);
            else % no comparions needed
                Yrand_opt = Yrand_all{1};
            end     
            
            [tmu, out_param] = meanMC_g(Yrand_opt, obj.in_param);
            out_param.ntot = out_param.ntot + nextra; % add counts of extra samples used
            out_param.optMethod = opt; % index of the optimal methods
            out_param.time = toc(tstart); % get running time n
        end
        
    end
    
    
    methods (Access = protected)
        
        function Y = genYrand_cv(obj, ridge)
            n = obj.cv_param.ncv * length(obj.cv_param.muX);
            YX = obj.cv_param.YXrand(n);% generate the matrix YX
            Y =YX(:,1); % get Y
            X = bsxfun(@minus,YX(:,2:end),mean(YX(:,2:end),1)); % get centered X
            if ridge == 0 % no regularization
                beta = bsxfun(@minus,X,mean(X,1))\Y; % get estimated coefficients
            else
                beta = (X'*X+eye(length(obj.cv_param.muX))*ridge)\X'*Y;
            end
            function fn = Yrand_cv(n) % construct new function
                YX_b = obj.cv_param.YXrand(n);
                Y_b = YX_b(:,1);
                X_b = YX_b(:,2:end);
                fn = Y_b - bsxfun(@minus,X_b,obj.cv_param.muX) * beta;
            end
            Y = @Yrand_cv;
        end
        
        function Y = genYrand_av(obj)
            function fn = Yrand_av(n) % construct new function
                YY = obj.av_param.YYrand(n);
                fn = (YY(:,1) + YY(:,2))/2; % average the antithetic pairs
            end
            Y = @Yrand_av;
        end
        
        
        function [Y_opt, index_opt] = selectYrand(obj, Yrand_all)
            num = length(Yrand_all); % number of methods to be compared
            time = zeros(1, num); % time to run each method
            variance = zeros(1, num); % estimated varaince for each method
            
            for i = 1:num
                fn = Yrand_all{i};
                tic, y = fn(obj.nc); time(i)=toc;
                variance(i) = var(y);
            end
            
            score = time .* variance; % score for each method,
            % roughly proportional to the total
            % time needs to estimate mu
            score_min = min(score); % lowest score
            index_opt  = (score == score_min);      
            Y_opt = Yrand_all{index_opt}; % return the method with the lowest score
        end
        
    end
    
    methods (Static, Access = protected)
        
        function [tmu,out_param]=meanMC_g(Yrand,in_param)
            tstart = tic;
            out_param = in_param;
            n1 = 2;
            Yrand(n1); %let it run once to load all the data. warm up the machine.
            nsofar = n1;
            
            ntry = 10;% try several samples to get the time
            ttry=time(@()Yrand(ntry));
            tpern = ttry/ntry; % calculate time per sample
            nsofar = nsofar+ntry; % update n so far
            out_param.exit = 0;
            if tpern<1e-7;%each sample use very very little time
                booster = 8;
                ttry2 = timeit(@()Yrand(ntry*booster));
                ntry = ntry*[1 booster];
                ttry = [ttry ttry2];% take eight times more samples to try
            elseif tpern>=1e-7 && tpern<1e-5 %each sample use very little time
                booster = 6;
                ttry2 = timeit(@()Yrand(ntry*booster));
                ntry = ntry*[1 booster];
                ttry = [ttry ttry2];% take six times more samples to try
            elseif tpern>=1e-5 && tpern<1e-3 %each sample use little time
                booster = 4;
                ttry2 = timeit(@()Yrand(ntry*booster));
                ntry = ntry*[1 booster];
                ttry = [ttry ttry2];% take four times more samples to try
            elseif  tpern>=1e-3 && tpern<1e-1 %each sample use moderate time
                booster = 2;
                ttry2 = timeit(@()Yrand(ntry*booster));
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
            sig0up = out_param.fudge.*sig0;% upper bound on the standard deviation
            alpha_sig = out_param.alpha/2;% the uncertainty for variance estimation
            alphai = (out_param.alpha-alpha_sig)/(2*(1-alpha_sig));
            %uncertainty to do iteration
            out_param.kurtmax = (out_param.nSig-3)/(out_param.nSig-1) ...
                + ((alpha_sig*out_param.nSig)/(1-alpha_sig))...
                *(1-1/out_param.fudge^2)^2;
            %the upper bound on the modified kurtosis
            tol1 = ncbinv(out_param.n1,alphai,out_param.kurtmax);
            %tolerance for initial estimation
            out_param.tol(1) = sig0up*tol1;
            %the width of initial confidence interval for the mean
            i=1;
            npcmax = 1e6;%constant to do iteration and mean calculation
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
                toltype = 'max';
                % error type, see the function 'tolfun' at Algoithms/+gail/ directory
                % for more info
                theta  = 0;% relative error tolerance case
                deltaplus = (gail.tolfun(out_param.abstol,out_param.reltol,...
                    theta,out_param.hmu(i) - out_param.tol(i),toltype)...
                    +gail.tolfun(out_param.abstol,out_param.reltol,...
                    theta,out_param.hmu(i) + out_param.tol(i),toltype))/2;
                % a combination of tolfun, which used to decide stopping time
                if deltaplus >= out_param.tol(i) % stopping criterion
                    deltaminus= (gail.tolfun(out_param.abstol,out_param.reltol,...
                        theta,out_param.hmu(i) - out_param.tol(i),toltype)...
                        -gail.tolfun(out_param.abstol,out_param.reltol,...
                        theta,out_param.hmu(i) + out_param.tol(i),toltype))/2;
                    % the other combination of tolfun, which adjust the hmu a bit
                    tmu = out_param.hmu(i)+deltaminus;
                    break;
                else
                    i=i+1;
                    deltat=0.7;
                    deltah=0.5;
                    delta=0;% constant to decide the next tolerance
                    out_param.tol(i) = max(min(deltaplus*deltat, ...
                        deltah*out_param.tol(i-1)),delta*out_param.tol(i-1));
                    %update the next tolerance
                    toloversig = out_param.tol(i)/sig0up;%next tolerance over sigma
                    alphai = (out_param.alpha-alpha_sig)/(1-alpha_sig)*2.^(-i);
                    %update the next uncertainty
                    out_param.n(i) = nchebe(toloversig,alphai,out_param.kurtmax);
                    %get the next sample size needed
                end
            end
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
        
        function tol1 = ncbinv(n1,alpha1,kurtmax)
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
            tol1 = min(NCheb_inv,NBE_inv);
            %take the min of Chebyshev and Berry Esseen tolerance
        end
        
        function meanMC_g_err(out_param)
            % Handles errors in meanMC_g and meanMC_g_param to give an exit with
            %  information.
            %            out_param.exit = 0   success
            %                             1   too many samples required
            
            if ~isfield(out_param,'exit'); return; end
            if out_param.exit==0; return; end
            switch out_param.exit
                case 1 % not enough samples to estimate the mean.
                    nexceed = out_param.n(out_param.tau);
                    warning('MATLAB:meanMC_g:maxreached',...
                        [' In order to achieve the guaranteed accuracy, at step '...
                        int2str(out_param.tau) ', tried to evaluate at ' int2str(nexceed) ...
                        ' samples, which is more than the remaining '...
                        int2str(out_param.nremain) ...
                        ' samples. We will use all the samples left to estimate the mean without guarantee.']);
                    return
            end
        end
        
    end
end
