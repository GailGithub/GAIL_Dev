function [Bpath,out_sample]=brownmotion(varargin)
% BROWNMOTION generates discrete Brownian motion paths
%
%   Bpath = BROWNMOTION(N) 
%   generates N paths of a standard Brownian motion on [0,1] with 8 time
%   steps.
%
%   Bpath = BROWNMOTION(N,T,d,disctype,randtype,s,check) 
%   generates N Brownian motion paths on [0,T] with d time steps.
%   Furthermore,
%       disctype = discretization type
%                  'timestep' = forward time step (default)
%                  'BB' =  Brownian bridge
%                  'KL' = Karhunen-Loeve
%       randtype = type of random vectors
%                  'IID' = independent and identically distributed (default)
%                  'Sobol' = Sobol points (a quasi-random sequence)
%              s = number of terms in BB or KL (default = 8)
%      needcheck = logical variable that determines whether parameters
%                  should be checked for validity (default=true)
%   The list of inputs can be truncated, and those missing are set to
%   default values.
%
%   Bpath = BROWNMOTION(N,'d',d,'randtype','Sobol',...)
%   generates N Brownian motion paths with parameters entered in arbitrary
%   order.
%
%   Bpath = BROWNMOTION(N,in_sample) 
%   generates N Brownian motion paths with optional parameters entered as
%   fields in the structure in_sample.  If a field is not specified, the
%   default value is used.

%% Parse and check the validity of input parameters
out_sample=bm_param(varargin{:});
N=out_sample.N;
T=out_sample.T;
d=out_sample.d;
disctype=out_sample.disctype;
s=out_sample.s;

%% Generate random vectors
switch disctype
    case 'timestep'
        dim=d; %number of time steps
    case 'BB'
        dim=s;
        d=s;
        out_sample.d=d;
    case 'KL'
        dim=s;
end

switch out_sample.randtype %generate Gaussian random vectors
    case 'IID'
        x=randn(N,dim);
    case 'Sobol'
        if ~isfield(out_sample,'Sobolstream')
            scsobol=qrandstream(scramble(sobolset(dim),...
                'MatousekAffineOwen')); %set up Sobol stream
            out_sample.Sobolstream=scsobol;
        end
        x=norminv(rand(out_sample.Sobolstream,N,dim));
            %get normal Sobol vectors
end

%% Turn random vectors into sample paths
delta=T/d;
tovT=(delta:delta:T);
out_sample.tnode=[0 tovT];

switch disctype
    case 'timestep' %forward time step
        Bpath=cumsum([zeros(N,1) x*sqrt(delta)],2);
    case 'BB' %Brownian Bridge
        log2s=log2(s);
        if floor(log2s)==log2s % do it fast
            saw=@(t) 1-abs(1-mod(t,2)); %triangular saw function
            Bpath=zeros(N,d+1);
            srange=[0 2.^(0:log2s)];
            slow=srange(1:log2s+1)+1;
            shi=srange(2:log2s+2);
            snum=diff(srange);
            factor=sqrt(T./[1 4*(shi(1:log2s))]);
            for k=1:log2s+1
                which=repmat(slow(k):shi(k),s/snum(k),1);
                which=which(:)';
                Bpath(:,2:d+1)=Bpath(:,2:d+1)+...
                    x(:,which).*repmat(saw(tovT*shi(k))*factor(k),N,1);
            end
        else %do it slowly
            hat=@(t) max(1-abs(t),0); %triangular hat function
            tcent=1-net(sobolset(1),s); %van der Corput sequence
            out_sample.tnode=[0 sort(tcent)]; %different nodes
            Bpath=zeros(N,d+1);       
            Bpath(:,2:d+1)=Bpath(:,2:d+1)+...
                x(:,1)*(tovT*sqrt(T));
            if s>1;
                tpowm=2.^(1+floor(log2(1:s-1)));
                for k=2:s;
                    Bpath(:,2:d+1)=Bpath(:,2:d+1)+...
                        x(:,k)*...
                        hat((tovT-tcent(k))*tpowm(k-1))*...
                        (sqrt(T/(tpowm(k-1)*2)));
                end
            end
        end
    case 'KL' %Karhunen-Loeve
        tovT=(delta:delta:T);
        Bpath=zeros(N,d+1);
        km1=pi*((1:s)-1/2);
        for k=1:s;
            Bpath(:,2:d+1)=Bpath(:,2:d+1)+ ...
                x(:,k)*(sin(km1(k)*tovT)*(sqrt(2*T)/km1(k)));
        end
end
end

function out_sample=bm_param(varargin)
% parse the input to the brownmotion.m function

%% Default parameter values
default.N = 1e4;% default number of paths
default.disctype = 'timestep';% default way to get Brownian motion
default.randtype = 'IID';% default type of random numbers
default.T  = 1;% default time horizon
default.d = 8;% default number of time steps
default.s = 8;% default number of bases
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
    in2=varargin{2};
    validvarargin=(isnumeric(in2) || isstruct(in2) ...
        || ischar(in2));
end

if ~validvarargin
    %only one valid input N, use all the default parameters
    out_sample.N=N;
    out_sample.T = default.T;
    out_sample.d = default.d;
    out_sample.disctype = default.disctype;
    out_sample.randtype = default.randtype;
    out_sample.s = default.s;
    out_sample.needcheck = default.needcheck;
else %parse the second input
    p = inputParser;
    addRequired(p,'N',@isnumeric);
    if isnumeric(in2)%if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'T',default.T,@isnumeric);
        addOptional(p,'d',default.d,@isnumeric);
        addOptional(p,'disctype',default.disctype,@isstr);
        addOptional(p,'randtype',default.randtype,@isstr);
        addOptional(p,'s',default.s,@isnumeric);
        addOptional(p,'needcheck',default.needcheck,@islogical);
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'T',default.T,@isnumeric);
        addParamValue(p,'d',default.d,@isnumeric);
        addParamValue(p,'disctype',default.disctype,@isstr);
        addParamValue(p,'randtype',default.randtype,@isstr);
        addParamValue(p,'s',default.s,@isnumeric);
        addParamValue(p,'needcheck',default.needcheck,@islogical);
    end
    parse(p,N,varargin{2:end})
    out_sample = p.Results;
end

%% Check parameter validity
if out_sample.needcheck
    if ~isposint(out_sample.N) % number of sample paths should be an integer
        warning(['The number of sample paths should be a positive integer, '...
            'using ' num2str(default.N)])
        out_sample.N = default.N;
    end

    if out_sample.T <= 0
        warning(['Time horizon should be greater than 0, ' ...
            'using ' num2str(default.T)])
        out_sample.T = default.T;
    end

    if ~ischar(out_sample.disctype)
        warning('Discretization type should be a character.')
        out_sample.disctype=[];
    end
    if any(strcmpi(out_sample.disctype,{'timestep','discrete'}))
        out_sample.disctype='timestep';
    elseif any(strcmpi(out_sample.disctype,{'BB','bridge'}))
        out_sample.disctype='BB';
    elseif any(strcmpi(out_sample.disctype,{'KL','Karhunen'}))
        out_sample.disctype='KL';
    else 
        warning(['Discretization type not recognized, using ' default.disctype])
        out_sample.disctype=default.disctype;
    end

    if any(strcmp(out_sample.disctype,{'timestep','KL'})) ...
            && ~isposint(out_sample.d) % number of time steps should be an integer
        warning(['The number of time steps should be a positive integer, '...
            'using ' int2str(default.d)])
        out_sample.d = default.d;
    end

    if any(strcmp(out_sample.disctype,{'BB','KL'})) ...
            && ~isposint(out_sample.s) % number of terms in expansion should be an integer
        warning(['The number of terms in the expansion should be a positive integer, '...
            'using ' int2str(default.s)])
        out_sample.s = default.s;
    end

    if ~ischar(out_sample.randtype)
        warning('Random vector type should be a character.')
        out_sample.randtype=[];
    end
    if any(strcmpi(out_sample.randtype,{'IID','random'}))
        out_sample.randtype='IID';
    elseif any(strcmpi(out_sample.randtype,{'Sobol','quasi-random'}))
        out_sample.randtype='Sobol';
    else 
        warning(['Discretization type not recognized, using ' default.randtype])
        out_sample.randtype=default.randtype;
    end
end 
end   
