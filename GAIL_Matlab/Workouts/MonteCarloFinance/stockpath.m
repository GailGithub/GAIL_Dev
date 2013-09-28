function [stock,out_stparam]=stockpath(varargin)
% STOCKPATH generates stock price paths for certain models
%
%   stock = STOCKPATH(N,in_stparam) 
%   generates N stock paths with the parameters listed above and others
%   provided as fields in the structure in_stparam.  If a field is not
%   specified, the default value is used.
%
%   The fields in in_stparam are
%       pathtype = type of stock price path
%                  'GBM' = geometric Brownian motion (default)       
%       T        = number of years to maturity (default = 1)
%       d        = number of time steps (default = 8)
%       S0       = initial stock price (default = 100)
%       r        = interest rate (default = 2%)
%       sig      = volatility (default = 50%)

out_stparam=stock_param(varargin{:});

N=out_stparam.N;
pathtype=out_stparam.pathtype;

T=out_stparam.T;
d=out_stparam.d;
time=(0:T/d:T);
S0=out_stparam.S0;
r=out_stparam.r;
sig=out_stparam.sig;

%% Generate stock paths
switch pathtype
    case 'GBM'
        [bmpath,out_stparam]=brownmotion(N,out_stparam);
        stock=S0*exp((r-sig.^2/2)*repmat(time,N,1)+sig*bmpath);
end
end

function out_stparam=stock_param(varargin)
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

%% Parse inputs
if isempty(varargin) %sample path number not input
    help brownmotion
    warning(['At least N must be specified. Now GAIL is using N = ' ...
        int2str(default.N)])
    N=default.N;
else
    N=varargin{1};
end

validvarargin=numel(varargin)>1;
if validvarargin
    in_stparam=varargin{2};
    validvarargin=isstruct(in_stparam);
end

if ~validvarargin
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
else %parse the structure in_stparam
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
    out_stparam.N=N;
end

%% Check parameter validity
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