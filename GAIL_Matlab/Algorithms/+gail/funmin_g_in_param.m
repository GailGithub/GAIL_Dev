% funmin_g_in_param: funmin_g's input parameter object
%
% Examples
%
% >> in_param = gail.funmin_g_in_param()
%    Warning: Function f must be a function handle. ***
%             f: @(x)exp(-100*(x-0.5).^2)
%             a: 0
%             b: 1
%        abstol: 1.0000e-06
%           nlo: 10
%           nhi: 1000
%          nmax: 10000000
%         ninit: 100
%
%
% >> f = @(x) x.^2; in_param = gail.funmin_g_in_param(f)
%    in_param = ***
%             f: @(x)x.^2
%             a: 0
%             b: 1
%        abstol: 1.0000e-06
%           nlo: 10
%           nhi: 1000
%          nmax: 10000000
%         ninit: 100
%
%
%  To get a struct:
%  >> in_param = gail.funmin_g_in_param( @(x) x.^2); out_param = in_param.toStruct()
%  out_param = 
% 
%             f: @(x)x.^2
%             a: 0
%             b: 1
%        abstol: 1.0000e-06
%           nlo: 10
%           nhi: 1000
%         ninit: 100
%          nmax: 10000000
%
%
% To get a structure with selected fields (and ignore properties that do not exist):
% >> in_param = gail.funmin_g_in_param( @(x) x.^2); out_param = in_param.toStruct({'f', 'abstol','c'})
%  out_param =
%
%         f: @(x)x.^2
%    abstol: 1.0000e-06
%
classdef funmin_g_in_param < gail.gail1D_in_param 
    %% data
    properties % public
        %TolX
    end % properties
     
    %% methods
    methods % public
        % constructor
        function out_param = funmin_g_in_param(varargin)
            % parse the input to a gail function
            in = cell(0);
            if nargin >= 1
                in = varargin{1};
            end
            out_param = out_param@gail.gail1D_in_param(in);
            out_param = out_param.set_input_field_names(...
              {'a','b','abstol',...%'TolX',
               'nlo','nhi','nmax'}...
            );
            % out_param.get_input_field_names();
            %% Default parameter values
            default = out_param.get_default();
            %default.TolX = 1e-3;
        
            %% parse inputs
            out_param = out_param.parse_inputs(default, varargin{:});
            
            %% validate inputs
            out_param = out_param.validate_inputs();
            
        end % constructor

        function out_param = validate_inputs(out_param)
            out_param = validate_inputs@gail.gail1D_in_param(out_param);
            
            %if (~isempty(out_param.TolX))
            %    if out_param.TolX < 0
            %        warning('GAIL:funmin_g_in_param:Xtolnonpos', ['X tolerance should be greater than or equal to 0.' ...
            %            ' Using default X tolerance ' num2str(default.TolX)]);
            %        out_param.TolX = default.TolX;
            %    end
            %end
        end
        
       function out_param = toStruct(out_param,varargin)
            l = {'f', 'a', 'b','abstol',...%'TolX',
                 'nlo','nhi','ninit','nmax','tau'};
            if length(varargin) > 0   
                l = varargin{1};
            end
            out_param = toStruct@gail.gail1D_in_param(out_param, l);
        end

    end % methods
    
    methods (Access = protected) % seen by subclasses
        
        % customize display order of properties (data) in an instance
        function propgrp = getPropertyGroups(~)
            proplist = {'f', 'a', 'b','abstol','nlo','nhi','nmax','ninit'}; % 'TolX','tau' 
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end % methods (protected)
    
    
end % classdef

%% other functions