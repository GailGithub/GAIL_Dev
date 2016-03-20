% funappx_g_in_param: GAIL funappx_g's input parameter object
%
% Examples
%
% >> in_param = gail.funappx_g_in_param()
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
%         exitflag: []
%             iter: []
%          npoints: []
%           errest: []
%
%
% >> f = @(x) x.^2; in_param = gail.funappx_g_in_param(f)
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
%         exitflag: []
%             iter: []
%          npoints: []
%           errest: []
%
%
% >> f = @(x) x.^2; in_param = gail.funappx_g_in_param(f,0,1,1e-6,10,1000,10000000,1000)
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
%         exitflag: []
%             iter: []
%          npoints: []
%           errest: []
%
%
% >> f = @(x) x.^2; in_param.a=0; in_param.b =1;  in_param = gail.funappx_g_in_param(f,in_param)
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
%         exitflag: []
%             iter: []
%          npoints: []
%           errest: []
%
%
%  To get a struct:
%  >> out_param = in_param.toStruct()
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
%     exitflag: []
%         iter: []
%      npoints: []
%       errest: []
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
classdef funappx_g_in_param < gail.gail_in_param & matlab.mixin.CustomDisplay
    %% data
    properties
        %% inputs
        a
        b
        abstol
        maxiter
        nhi
        nlo
        nmax
        nstar
        ninit
        
        % optional
        memorytest
        output_x
        
        %% outputs
        exitflag
        errest
        npoints
        iter
        
    end % properties
    
    properties (GetAccess = private, SetAccess = private)  % seen by subclasses
        input_field_names = {'a','b','abstol','nlo','nhi','nmax','maxiter','memorytest','output_x'};
        output_field_names = {'ninit','exitflag','iter','npoints','errest'};
    end
    
    %% methods
    methods % public
        
        % constructor
        function out_param = funappx_g_in_param(varargin)
            %% parse the input to a gail function
            if nargin >= 1
                in = varargin{1};
            else
                in = cell(0);
            end
            out_param = out_param@gail.gail_in_param(in);
            %% Default parameter values
            default.a = 0;
            default.b = 1;
            default.abstol = 1e-6;
            default.nlo = 10;
            default.nhi = 1000;
            default.nmax = 1e7;
            default.maxiter = 1000;
            default.memorytest = 0;
            default.output_x = 0;
            
            %% parse inputs
            out_param = out_param.parse_inputs(default, varargin{:});
            
            %% validate inputs
            out_param = validate_inputs(out_param);
        end % constructor
        
        function out_param = parse_inputs(out_param,default,varargin)
            MATLABVERSION = gail.matlab_version;
            if MATLABVERSION >= 8.3
                f_addParamVal = @addParameter;
            else
                f_addParamVal = @addParamValue;
            end;
            
            validvarargin=numel(varargin)>1;
            if validvarargin
                in2=varargin{2};
                validvarargin=(isnumeric(in2) || isstruct(in2) ...
                    || ischar(in2));
            end
            
            if ~validvarargin
                %if only one input f, use all the default parameters
                for field_index = 1:length(out_param.input_field_names)
                    field = out_param.input_field_names{field_index};
                    out_param.(field) = default.(field);
                end
            else
                p = inputParser;
                addRequired(p,'f',@gail.isfcn);
                if isnumeric(in2)%if there are multiple inputs with
                    %only numeric, they should be put in order.
                    for field_index = 1:length(out_param.input_field_names)
                        field = out_param.input_field_names{field_index};
                        addOptional(p,field,default.(field),@isnumeric);
                    end
                else
                    if isstruct(in2) %parse input structure
                        p.StructExpand = true;
                        p.KeepUnmatched = true;
                    end
                    for field_index = 1:length(out_param.input_field_names)
                        field = out_param.input_field_names{field_index};
                        f_addParamVal(p,field,default.(field),@isnumeric);
                    end
                end
                parse(p,out_param.f,varargin{2:end})
                results = p.Results;
                field_list = fields(results);
                for field_index = 1:length(field_list)
                    field = field_list{field_index};
                    out_param.(field) = results.(field);
                end
            end;
        end % function parse_input
        
        function out_param = validate_inputs(out_param)
            % let end point of interval not be infinity
            if (out_param.a == inf||out_param.a == -inf)
                warning('GAIL:funappx_g_in_param:aisinf',['a cannot be infinity. '...
                    'Use default a = ' num2str(default.a)])
                out_param.a = default.a;
            end;
            if (out_param.b == inf||out_param.b == -inf)
                warning(['GAIL:funappx_g_in_param:bisinf','b cannot be infinity. '...
                    'Use default b = ' num2str(default.b)])
                out_param.b = default.b;
            end;
            
            if (out_param.b < out_param.a)
                warning('GAIL:funappx_g_in_param:blea',['b cannot be smaller than a;'...
                    ' exchange these two. '])
                tmp = out_param.b;
                out_param.b = out_param.a;
                out_param.a = tmp;
            elseif(out_param.b == out_param.a)
                warning('GAIL:funappx_g_in_param:beqa',['b cannot equal a. '...
                    'Use b = ' num2str(out_param.a+1)])
                out_param.b = out_param.a+1;
            end;
            
            % let error tolerance greater than 0
            if (out_param.abstol <= 0 )
                warning('GAIL:funappx_g_in_param:tolneg', ['Error tolerance should be greater'...
                    ' than 0. Using default error tolerance ' num2str(default.abstol)])
                out_param.abstol = default.abstol;
            end
            % let cost budget be a positive integer
            if (~gail.isposint(out_param.nmax))
                if gail.isposintive(out_param.nmax)
                    warning('GAIL:funappx_g_in_param:budgetnotint',['Cost budget should be '...
                        'a positive integer. Using cost budget '...
                        , num2str(ceil(out_param.nmax))])
                    out_param.nmax = ceil(out_param.nmax);
                else
                    warning('GAIL:funappx_g_in_param:budgetisneg',['Cost budget should be '...
                        'a positive integer. Using default cost budget '...
                        int2str(default.nmax)])
                    out_param.nmax = default.nmax;
                end;
            end
            
            if (~gail.isposint(out_param.nlo))
                if gail.isposge3(out_param.nlo)
                    warning('GAIL:funappxglobal_g:lowinitnotint',['Lower bound of '...
                        'initial number of points should be a positive integer.' ...
                        ' Using ', num2str(ceil(out_param.nlo)) ' as nlo '])
                    out_param.nlo = ceil(out_param.nlo);
                else
                    warning('GAIL:funappxglobal_g:lowinitlt3',[' Lower bound of '...
                        'initial number of points should be a positive integer greater'...
                        ' than 3. Using 3 as nlo'])
                    out_param.nlo = 3;
                end
                warning('GAIL:funappxglobal_g:lowinitnotint',['Lower bound of '...
                    'initial nstar should be a positive integer.' ...
                    ' Using ', num2str(ceil(out_param.nlo)) ' as nlo '])
                out_param.nlo = ceil(out_param.nlo);
            end
            if (~gail.isposint(out_param.nhi))
                if gail.isposge3(out_param.nhi)
                    warning('GAIL:funappxglobal_g:hiinitnotint',['Upper bound of '...
                        'initial number of points should be a positive integer.' ...
                        ' Using ', num2str(ceil(out_param.nhi)) ' as nhi' ])
                    out_param.nhi = ceil(out_param.nhi);
                else
                    warning('GAIL:funappxglobal_g:hiinitlt3',[' Upper bound of '...
                        'points should be a positive integer greater than 3. Using '...
                        'default number of points ' int2str(default.nhi) ' as nhi' ])
                    out_param.nhi = default.nhi;
                end
                warning('GAIL:funappxglobal_g:hiinitnotint',['Upper bound of '...
                    'initial nstar should be a positive integer.' ...
                    ' Using ', num2str(ceil(out_param.nhi)) ' as nhi' ])
                out_param.nhi = ceil(out_param.nhi);
            end
            
            if (out_param.nlo > out_param.nhi)
                warning('GAIL:funappx_g_in_param:logrhi', ['Lower bound of initial number of'...
                    ' points is larger than upper bound of initial number of '...
                    'points; Use nhi as nlo'])
                out_param.nhi = out_param.nlo;
            end;
            
            h = out_param.b - out_param.a;
            out_param.ninit = ceil(out_param.nhi*(out_param.nlo/out_param.nhi)...
                ^(1/(1+h)));
            
            if (~gail.isposint(out_param.maxiter))
                if gail.ispositive(out_param.maxiter)
                    warning('GAIL:funappx_g_in_param:maxiternotint',['Max number of '...
                        'iterations should be a positive integer. Using max number '...
                        'of iterations as  ', num2str(ceil(out_param.maxiter))])
                    out_param.nmax = ceil(out_param.nmax);
                else
                    warning('GAIL:funappx_g_in_param:budgetisneg',['Max number of iterations'...
                        ' should be a positive integer. Using max number of '...
                        'iterations as ' int2str(default.maxiter)])
                    out_param.nmax = default.nmax;
                end;
            end
            if (out_param.memorytest~=1&&out_param.memorytest~=0)
                warning('GAIL:funappx_g_in_param:memorytest', ['Input of memorytest'...
                    ' can only be 1 or 0; use default value 0'])
                out_param.memorytest = 0;
            end;
            if (out_param.output_x~=1&&out_param.output_x~=0)
                warning('GAIL:funappx_g_in_param:output_x', ['Input of output_x'...
                    ' can only be 1 or 0; use default value 0'])
                out_param.output_x = 0;
            end;
            
        end % function validate_inputs
        
        
        function out_struct = toStruct(out_param,varargin)
            field_list = ...%union({'f'}, union(out_param.input_field_names, out_param.output_field_names,'stable'),'stable');
                {'f','a','b','abstol','nlo','nhi','ninit','nmax','maxiter',...
                'exitflag','iter','npoints','errest'};
            if nargin > 1
                field_list = varargin{1};
            end
            for field_index = 1:length(field_list)
                field = field_list{field_index};
                if isprop(out_param,field)
                    out_struct.(field) = out_param.(field);
                end
            end
        end
        
        function names = get_input_field_names(out_param)
            names = out_param.input_field_names;
        end
        
        function names = get_output_field_names(out_param)
            names = out_param.output_field_names;
        end
    end % methods
    
    
    methods (Access = protected)  % seen by subclasses
        
        % customize display order of properties (data) in an instance
        function propgrp = getPropertyGroups(~)
            proplist = {'f', 'a', 'b','abstol','nlo','nhi','ninit','nmax','maxiter',...
                'exitflag','iter','npoints','errest','x'};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end % methods (protected)
    
    
end % classdef

%% other functions