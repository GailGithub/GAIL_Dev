function [pay,out_payparam,out_stparam]=payoff_XS(varargin)
% PAYOFF generates option payoffs for certain options
%
%   pay = PAYOFF(N,in_payparam,in_stparam) 
%   generates N payoffs of a European call option with the default
%   parameters
%       stock    = function handle for stock price paths 
%                  (default = GBM)
%       paytype  = type of option
%                  'eurocall' = European call (default)
%                  'europut' = European put
%                  'ameancall' = European call
%                  'ameanput' = European put
%          
%       T        = number of years to maturity (default = 1)
%       K        = strike price (default = 100)
%       r        = interest rate (default = 2%)
%       barrier  = option barrier (default = 120)
%
%   pay = PAYOFF(N,in_payparam,in_stparam) generates N payoffs of an option
%   with the fields in in_payparam determining the payoff and the fields in
%   in_stparam determining the stock paths.

[out_payparam,out_stparam]=pay_param(varargin{:});
N=out_stparam.N;
T=out_stparam.T;
d=out_stparam.d;
r=out_stparam.r;
paytype=out_payparam.paytype;
K=out_payparam.K;
barrier=out_payparam.barrier;

%% Generate random vectors
Ndmax=1e7;
Nmax=floor(Ndmax/d);
pc=ceil(N/Nmax);
Nrest=N-(pc-1)*Nmax;
Npc=[repmat(Nmax,pc-1,1); Nrest];
Nvec=[0; cumsum(Npc)];
pay=zeros(N,1);
temp_stparam=out_stparam;
for j=1:pc
    n=Npc(j);
    stock=stockpath(n,temp_stparam);
    switch paytype
        case 'eurocall'
            pay(Nvec(j)+1:Nvec(j+1))=max(stock(:,d+1)-K,0)*exp(-r*T);
        case 'europut'
            pay(Nvec(j)+1:Nvec(j+1))=max(K-stock(:,d+1),0)*exp(-r*T);
        case {'ameancall','ameanput','gmeancall','gmeanput'}
            if strcmp(paytype(1),'a') %arithmetic mean
                meanstock=mean(stock(:,2:d+1),2);
            else %geometric mean
                meanstock=(prod(stock(:,2:d+1),2)).^(1/d);
            end
            %length(Nvec(j)+1:Nvec(j+1)), length(meanstock), keyboard
            if strcmp(paytype(end-3:end),'call')
                pay(Nvec(j)+1:Nvec(j+1))=max(meanstock-K,0)*exp(-r*T);
            else
                pay(Nvec(j)+1:Nvec(j+1))=max(K-meanstock,0)*exp(-r*T);
            end
        %Barrier option
        case {'upincall','upoutcall','downincall','downoutcall','upinput','upoutput','downinput','downoutput'}
            if strcmp(paytype(1:2),'up') %up barrier
                in=any(stock>barrier);
                if strcmp(paytype(3:4),'in') %up and in barrier
                    if strcmp(paytype(end-3:end),'call') %up and in eurocall
                        pay(Nvec(j)+1:Nvec(j+1))=max(stock(:,d+1)-K,0)*exp(-r*T)*in;
                    else %up and in europut
                        pay(Nvec(j)+1:Nvec(j+1))=max(K-stock(:,d+1),0)*exp(-r*T)*in;
                    end
                else %up and out barrier
                    if strcmp(paytype(end-3:end),'call') %up and out eurocall
                        pay(Nvec(j)+1:Nvec(j+1))=max(stock(:,d+1)-K,0)*exp(-r*T)*~in;
                    else %up and out europut
                        pay(Nvec(j)+1:Nvec(j+1))=max(K-stock(:,d+1),0)*exp(-r*T)*~in;
                    end
                end
            else %down barrier
                in=any(stock<barrier);
                if strcmp(paytype(5:6),'in') %down and in barrier
                    if strcmp(paytype(end-3:end),'call') %down and in eurocall
                        pay(Nvec(j)+1:Nvec(j+1))=max(stock(:,d+1)-K,0)*exp(-r*T)*in;
                    else %down and in europut
                        pay(Nvec(j)+1:Nvec(j+1))=max(K-stock(:,d+1),0)*exp(-r*T)*in;
                    end
                else %down and out barrier
                    if strcmp(paytype(end-3:end),'call') %down and out eurocall
                        pay(Nvec(j)+1:Nvec(j+1))=max(stock(:,d+1)-K,0)*exp(-r*T)*~in;
                    else %down and out europut
                        pay(Nvec(j)+1:Nvec(j+1))=max(K-stock(:,d+1),0)*exp(-r*T)*~in;
                    end
                end
            end
        case 'lookcall'
            pay(Nvec(j)+1:Nvec(j+1))=(stock(:,d+1)-min(stock))*exp(-r*T);
        case 'lookput'
            pay(Nvec(j)+1:Nvec(j+1))=(max(stock)-stock(:,d+1))*exp(-r*T);
    end
end
end

function [out_payparam,out_stparam]=pay_param(varargin)
% parse the input to the brownmotion.m function

%% Default parameter values
default.N = 1e4;% default number of paths
default.T  = 1;% default time horizon
default.d  = 8;% default number of time steps
default.S0 = 100;% default initial stock price
default.r = 0.02;% default interest rate
default.sig = 0.5;% default volatility
default.pathtype='GBM';% default path type
default.randtype='IID';% default random numbers
default.disctype='timestep';% default random numbers
default.needcheck = true;% default of whether to check parameter validity
default.paytype = 'eurocall';% default way to get Brownian motion
default.K = 100;% default strike price
default.barrier = 120;% default barrier value

%% Parse inputs
if isempty(varargin) %sample path number not input
    help brownmotion
    warning(['At least N must be specified. Now GAIL is using N = ' ...
        int2str(default.N)])
    N=default.N;
else
    N=varargin{1};
end

validpaystr=numel(varargin)>1;
if validpaystr
    in_payparam=varargin{2};
    validpaystr=isstruct(in_payparam);
end

validstockstr=validpaystr&&(numel(varargin)>2);
if validstockstr
    in_stparam=varargin{3};
    validstockstr=isstruct(in_stparam);
end

if ~validpaystr
    %only one valid input N, use all the default parameters
    out_payparam.paytype = default.paytype;
    out_payparam.K = default.K;
    out_payparam.needcheck=default.needcheck;
    out_payparam.barrier=default.barrier;
else %parse the second input
    p = inputParser;
    addRequired(p,'N',@isnumeric);
    p.StructExpand = true;
    p.KeepUnmatched = true;
    addParamValue(p,'paytype',default.paytype,@isstr);
    addParamValue(p,'K',default.K,@isnumeric);
    addParamValue(p,'needcheck',default.needcheck,@islogical);
    addParamValue(p,'barrier',default.barrier,@isnumeric);
    parse(p,N,in_payparam)
    out_payparam = p.Results;
end

if ~validstockstr
    %only one valid input N, use all the default parameters
    out_stparam.N=N;
    out_stparam.T = default.T;
    out_stparam.d = default.d;
    out_stparam.S0 = default.S0;
    out_stparam.r = default.r;
    out_stparam.sig=default.sig;
    out_stparam.pathtype=default.randtype;
    out_stparam.randtype=default.randtype;
    out_stparam.disctype=default.disctype;
    out_stparam.needcheck=default.needcheck;
else %parse in_stparam
    p = inputParser;
    addRequired(p,'N',@isnumeric);
    p.StructExpand = true;
    p.KeepUnmatched = true;
    addParamValue(p,'T',default.T,@isnumeric);
    addParamValue(p,'d',default.d,@isnumeric);
    addParamValue(p,'S0',default.S0,@isnumeric);
    addParamValue(p,'r',default.r,@isnumeric);
    addParamValue(p,'sig',default.sig,@isnumeric);
    addParamValue(p,'pathtype',default.pathtype,@isstr);
    addParamValue(p,'disctype',default.disctype,@isstr);
    addParamValue(p,'randtype',default.randtype,@isstr);
    addParamValue(p,'needcheck',default.needcheck,@islogical);
    parse(p,N,in_stparam)
    out_stparam = p.Results;
end

%% Check parameter validity
if out_payparam.needcheck
    if ~ischar(out_payparam.paytype)
        warning('Option type should be a character.')
        out_payparam.paytype=[];
    end
    if ~any(strcmpi(out_payparam.paytype,...
            {'eurocall','europut',...
            'ameancall','ameanput'...
            'gmeancall','gmeanput'...
            'upincall','upoutcall','downincall','downoutcall'...
            'upinput','upoutput','downinput','downoutput'...
            'lookcall','lookput'}))
        warning(['Payoff type not recognized, using ' default.paytype])
        out_payparam.paytype='eurocall';
    end

    if out_payparam.K <= 0
        warning(['Strike price should be greater than 0, ' ...
            'using ' num2str(default.K)])
        out_payparam.K = default.K;
    end
    
    if out_payparam.barrier <= 0
        warning(['Barrier value should be greater than 0, ' ...
            'using ' num2str(default.barrier)])
        out_payparam.barrier = default.barrier;
    end
    
    out_payparam.needcheck=false;
end


if out_stparam.needcheck
    if ~isposint(out_stparam.N) % number of sample paths should be an integer
        warning(['The number of sample paths should be a positive integer, '...
            'using ' num2str(default.N)])
        out_stparam.N = default.N;
    end

    if out_stparam.T <= 0
        warning(['Time horizon should be greater than 0, ' ...
            'using ' num2str(default.T)])
        out_stparam.T = default.T;
    end

    if ~isposint(out_stparam.d)
        warning(['Number of time steps must be a positive integer , ' ...
            'using ' num2str(default.d)])
        out_stparam.d = default.d;
    end

    if out_stparam.S0 < 0
        warning(['Initial stock price should be no less than 0, ' ...
            'using ' num2str(default.S0)])
        out_stparam.S0 = default.S0;
    end

    if out_stparam.r < 0
        warning(['Interest rate should be no less than 0, ' ...
            'using ' num2str(default.r)])
        out_stparam.r = default.r;
    end

    if out_stparam.sig < 0
        warning(['Volatility should be greater than 0, ' ...
            'using ' num2str(default.sig)])
        out_stparam.sig = default.sig;
    end

    if ~ischar(out_stparam.pathtype)
        warning('Stock path type should be a character.')
        out_stparam.type=[];
    end
    if ~any(strcmpi(out_stparam.pathtype,...
            {'GBM'}))
        warning(['Stock path type not recognized, using ' default.pathtype])
        out_stparam.pathtype=default.pathtype;
    end

    if ~ischar(out_stparam.disctype)
        warning('Discretization type should be a character.')
        out_stparam.disctype=[];
    end
    if any(strcmpi(out_stparam.disctype,{'timestep','discrete'}))
        out_stparam.disctype='timestep';
    elseif any(strcmpi(out_stparam.disctype,{'BB','bridge'}))
        out_stparam.disctype='BB';
    elseif any(strcmpi(out_stparam.disctype,{'KL','Karhunen'}))
        out_stparam.disctype='KL';
    else 
        warning(['Discretization type not recognized, using ' default.disctype])
        out_stparam.disctype=default.disctype;
    end

    if any(strcmp(out_stparam.disctype,{'timestep','KL'})) ...
            && ~isposint(out_stparam.d) % number of time steps should be an integer
        warning(['The number of time steps should be a positive integer, '...
            'using ' int2str(default.d)])
        out_stparam.d = default.d;
    end

    if any(strcmp(out_stparam.disctype,{'BB','KL'})) ...
            && ~isposint(out_stparam.s) % number of terms in expansion should be an integer
        warning(['The number of terms in the expansion should be a positive integer, '...
            'using ' int2str(default.s)])
        out_stparam.s = default.s;
    end

    if ~ischar(out_stparam.randtype)
        warning('Random vector type should be a character.')
        out_stparam.randtype=[];
    end
    if any(strcmpi(out_stparam.randtype,{'IID','random'}))
        out_stparam.randtype='IID';
    elseif any(strcmpi(out_stparam.randtype,{'Sobol','quasi-random'}))
        out_stparam.randtype='Sobol';
    else 
        warning(['Discretization type not recognized, using ' default.randtype])
        out_stparam.randtype=default.randtype;
    end
    
    out_stparam.needcheck=false;

end

end   