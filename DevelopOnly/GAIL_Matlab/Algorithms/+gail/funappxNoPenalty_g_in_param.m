% funappxNoPenalty_g_in_param: funappxNoPenalty_g's input parameter object
%
% Examples
%
% >> in_param = gail.funappxNoPenalty_g_in_param()
%     Warning: Function f must be a function handle.***
%     in_param = ***
%                f: @(x)exp(-100*(x-0.5).^2)
%                a: 0
%                b: 1
%           abstol: 1.0000e-06
%              nlo: 10
%              nhi: 1000
%            ninit: 100
%             nmax: 10000000
%          maxiter: 1000
%
%
% >> f = @(x) x.^2; in_param = gail.funappxNoPenalty_g_in_param(f)
%    in_param = ***
%                f: @(x)x.^2
%                a: 0
%                b: 1
%           abstol: 1.0000e-06
%              nlo: 10
%              nhi: 1000
%            ninit: 100
%             nmax: 10000000
%          maxiter: 1000
%
%
% >> f = @(x) x.^2; in_param = gail.funappxNoPenalty_g_in_param(f,0,1,1e-6,10,1000,10000000,1000)
%    in_param = ***
%                f: @(x)x.^2
%                a: 0
%                b: 1
%           abstol: 1.0000e-06
%              nlo: 10
%              nhi: 1000
%            ninit: 100
%             nmax: 10000000
%          maxiter: 1000
%
%
% >> f = @(x) x.^2; in_param.a=0; in_param.b =1;  in_param = gail.funappxNoPenalty_g_in_param(f,in_param)
%    in_param = ***
%                f: @(x)x.^2
%                a: 0
%                b: 1
%           abstol: 1.0000e-06
%              nlo: 10
%              nhi: 1000
%            ninit: 100
%             nmax: 10000000
%          maxiter: 1000
%
%
%  To get a struct:
%  >> in_param = gail.funappxNoPenalty_g_in_param(@(x) x.^2); out_param = in_param.toStruct()
%   out_param =
%
%            f: @(x)x.^2
%            a: 0
%            b: 1
%       abstol: 1.0000e-06
%          nlo: 10
%          nhi: 1000
%        ninit: 100
%         nmax: 10000000
%      maxiter: 1000
%
%
% To get a structure with selected fields (and ignore properties that do not exist):
% >> in_param = gail.funappxNoPenalty_g_in_param(@(x) x.^2); out_param = in_param.toStruct({'f','abstol','c'})
%  out_param =
%
%     f: @(x)x.^2
%     abstol: 1.0000e-06
%
classdef funappxNoPenalty_g_in_param < gail.gail1D_in_param
    %% data
    properties % public
    end % properties
     
    %% methods
    methods % public
        % constructor
        function out_param = funappxNoPenalty_g_in_param(varargin)
            % parse the input to a gail function
            in = cell(0);
            if nargin >= 1
                in = varargin{1};
            end
            out_param = out_param@gail.gail1D_in_param(in);

            %% Default parameter values
            default = out_param.get_default();
        
            %% parse inputs
            out_param = out_param.parse_inputs(default, varargin{:});
            
            %% validate inputs
            out_param = out_param.validate_inputs();
            
        end % constructor

        
        function out_param = validate_inputs(out_param)
            out_param = validate_inputs@gail.gail1D_in_param(out_param);
            
            if (~isempty(out_param.abstol))
                if (out_param.abstol <= 0 )
                    warning('GAIL:funappxNoPenalty_g_in_param:tolneg', ['Error tolerance should be greater'...
                        ' than 0. Using default error tolerance ' num2str(default.abstol)])
                    out_param.abstol = default.abstol;
                end
            end
        end

    end % methods
    
    methods (Access = protected) % seen by subclasses
       
    end % methods (protected)

end % classdef

%% other functions