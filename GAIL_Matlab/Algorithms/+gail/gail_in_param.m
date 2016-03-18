% Prototype of GAIL input parameter object
%
% Examples
%
% >> in_param = gail_in_param()
%
%    in_param =
%      gail_in_param with properties:
%              f: @(x)exp(-100*(x-0.5).^2)
%              a: 0
%              b: 1
%         abstol: 1.0000e-06
%        maxiter: 1000
%            nhi: 1000
%            nlo: 10
%           nmax: 10000000
%          nstar: []
%          ninit: 100
%           exit: []
%           iter: []
%        npoints: []
%         errest: []
%     memorytest: 0
%       output_x: 0
%
%
% >> f = @(x) x.^2; in_param = gail_in_param(f)
%
%    in_param =
%      gail_in_param with properties:
%              f: @(x)x.^2
%              a: 0
%              b: 1
%         abstol: 1.0000e-06
%        maxiter: 1000
%            nhi: 1000
%            nlo: 10
%           nmax: 10000000
%          nstar: []
%          ninit: 100
%           exit: []
%           iter: []
%        npoints: []
%         errest: []
%     memorytest: 0
%       output_x: 0
%
classdef gail_in_param
    %% data
    properties
        % input function handle
        f    
        
    end
    
    %% methods
    methods
        % constructor
        function out_param = gail_in_param(varargin)
            
            % parse the input to a gail function
           
            MATLABVERSION = gail.matlab_version;
            if MATLABVERSION >= 8.3
                f_addParamVal = @addParameter;
            else
                f_addParamVal = @addParamValue;
            end;
            
            if isempty(varargin)
                warning('GAIL:gail_in_param:nofunction',['Function f must be specified. '...
                    'Now GAIL is using f(x)=exp(-100*(x-0.5)^2) and unit interval '...
                    '[0,1].'])
                help gail_in_param
                f = @(x) exp(-100*(x-0.5).^2);
                out_param.f = f;
            else
                if gail.isfcn(varargin{1})
                    f = varargin{1};
                    out_param.f = f;
                else
                    warning('GAIL:gail_in_param:notfunction',['Function f must be a '...
                        'function handle. Now GAIL is using f(x)=exp(-100*(x-0.5)^2).'])
                    f = @(x) exp(-100*(x-0.5).^2);
                    out_param.f = f;
                end
            end; 
        end % constructor
    end % methods
    
end % classdef

%% other functions