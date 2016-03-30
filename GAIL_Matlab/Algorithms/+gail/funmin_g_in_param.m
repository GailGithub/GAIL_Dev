% Prototype of GAIL input parameter object
%
% Examples
%
% >> in_param = gail.funmin_g_in_param()
%    Warning: Function f must be a function handle. ***
%             f: @(x)exp(-100*(x-0.5).^2)
%             a: 0
%             b: 1
%        abstol: 1.0000e-06
%          TolX: 1.0000e-03
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
%          TolX: 1.0000e-03
%           nlo: 10
%           nhi: 1000
%          nmax: 10000000
%         ninit: 100
%
%
%  To get a struct:
%  >> out_param = in_param.toStruct()
%  out_param = 
% 
%             f: @(x)x.^2
%             a: 0
%             b: 1
%        abstol: 1.0000e-06
%          TolX: 1.0000e-03
%           nlo: 10
%           nhi: 1000
%         ninit: 100
%          nmax: 10000000
%
%
% To get a structure with selected fields (and ignore properties that do not exist):
% >> out_param = in_param.toStruct({'f', 'a', 'b','c'})
%  out_param =
%
%     f: @(x)x.^2
%     a: 0
%     b: 1
%
classdef funmin_g_in_param < gail.gail1D_in_param 
    %% data
    properties % public
        TolX
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
              {'a','b','abstol','TolX','nlo','nhi','nmax'}...
            );
           % out_param.get_input_field_names();
            %% Default parameter values
            default = out_param.get_default();
            default.TolX = 1e-3;
        
            %% parse inputs
            out_param = out_param.parse_inputs(default, varargin{:});
            
            %% validate inputs
            
            out_param = out_param.validate_inputs();
            
        end % constructor
         
        function out_param = toStruct(out_param,varargin)
            out_param = toStruct@gail.gail1D_in_param(out_param,{'f', 'a', 'b','abstol','TolX','nlo','nhi','ninit','nmax','tau'});
                %'exitflag', 'npoints', 'errest','volumeX','tauchange','intervals'});
        end
        
        function out_param = validate_inputs(out_param)
            out_param = validate_inputs@gail.gail1D_in_param(out_param);
            
%             % Check whether the length tolerance is nonnegative
%             if out_param.abstol < 0
%                 warning('GAIL:funmin_g_in_param:abstolnonpos', ['Error tolerance should be greater than or equal to 0.' ...
%                     ' Using default error tolerance ', num2str(default.abstol)])
%                 out_param.abstol = default.abstol;
%             end
            
            if out_param.TolX < 0
                warning('GAIL:funmin_g_in_param:Xtolnonpos', ['X tolerance should be greater than or equal to 0.' ...
                    ' Using default X tolerance ' num2str(default.TolX)]);
                out_param.TolX = default.TolX;
            end
            

        end

    end % methods
    
    methods (Access = protected) % seen by subclasses
        
        % customize display order of properties (data) in an instance
        function propgrp = getPropertyGroups(~)
            proplist = {'f', 'a', 'b','abstol','TolX','nlo','nhi','nmax','ninit','tau'};
              %  'exitflag','npoints','errest','volumeX', 'tauchange', 'intervals'};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end % methods (protected)
    
    
end % classdef

%% other functions