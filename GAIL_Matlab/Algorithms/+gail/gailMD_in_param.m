% gailMD_in_param: GAIL input parameter object for multidimensional algorithms
%
% Examples
%
% >> in_param = gail.gailMD_in_param()
% Warning: Function f must be a function handle. Now GAIL is using f(x)=exp(-100*(x-0.5)^2). 
% > In gail_in_param>gail_in_param.gail_in_param at ***
%   In gailMD_in_param>gailMD_in_param.gailMD_in_param at ***
% 
% in_param = 
% 
%   gailMD_in_param with properties:
%                f: @(x)exp(-100*(x-0.5).^2)
%          measure: 'uniform'
%           abstol: 0.0100
%           reltol: 0.1000
%            fudge: 1.2000
%              dim: 1
%         hyperbox: [2x1 double]
%
%
%
% >> f = @(x) x.^2; in_param = gail.gailMD_in_param(f)
%     in_param = 
% 
%       gailMD_in_param with properties:
% 
%                f: @(x)x.^2
%          measure: 'uniform'
%           abstol: 0.0100
%           reltol: 0.1000
%            fudge: 1.2000
%              dim: 1
%         hyperbox: [2x1 double]
%
%
%
%
%  To get a struct:
%  >> in_param = gail.gailMD_in_param(@(x) x.^2); out_param = in_param.toStruct()
%   out_param =
%            f: @(x)x.^2
%      measure: 'uniform'
%       abstol: 0.0100
%       reltol: 0.1000
%        fudge: 1.2000
%          dim: 1
%     hyperbox: [2x1 double]
%
%
% To get a structure with selected fields (and ignore properties that do not exist):
% >> out_param = in_param.toStruct({'f', 'measure', 'hyperbox','nonexistent'})
%  out_param =
%            f: @(x)x.^2
%      measure: 'uniform'
%     hyperbox: [2x1 double]
%
classdef gailMD_in_param < gail.gail_in_param & matlab.mixin.CustomDisplay
    %% data
    properties
        %% inputs
        measure
        hyperbox
        dim
        abstol
        reltol
        fudge
    end % properties
    
    properties (GetAccess = protected, SetAccess = protected)  % seen by subclasses
        input_field_names = {'measure','abstol','reltol','fudge','dim','hyperbox'}; %order of parsing
        default_values
    end
    
    %% methods
    methods % public

        % constructor
        function out_param = gailMD_in_param(varargin)
            %% parse the input to a gail function
            if nargin >= 1
                in = varargin{1};
            else
                in = cell(0);
            end
            out_param = out_param@gail.gail_in_param(in);
            
            %% Default parameter values
            default.measure = 'uniform';% default measure
            default.dim = 1;% default dimension
            default.hyperbox = [zeros(1,default.dim);ones(1,default.dim)];% default hyperbox
            default.abstol  = 1e-2;% default absolute error tolerance
            default.reltol  = 1e-1;% default absolute error tolerance
            default.fudge = 1.2; % default variance inflation factor
   
            out_param.default_values = default;
            %% parse inputs
            out_param = out_param.parse_inputs(default, varargin{:});
            
            %% validate inputs
            out_param = out_param.validate_inputs();
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
                in2 = varargin{2};
                validvarargin = (isnumeric(in2) || isstruct(in2) ...
                    || ischar(in2));
            end
            
            field_names = out_param.input_field_names; %out_param.input_field_names; %;
            %out_param.get_input_field_names();
            if ~validvarargin
                %if only one input f, use all the default parameters
                for field_index = 1:length(field_names)
                    field = field_names{field_index};
                    out_param.(field) = default.(field);
                end
            else
                p = inputParser;
                addRequired(p,'f',@gail.isfcn);
                if isnumeric(in2)%if there are multiple inputs with
                    %only numeric, they should be put in order.
                    p.StructExpand = false;
                    for field_index = 1:length(field_names)
                        field = field_names{field_index};
                        addOptional(p,field,default.(field),@isnumeric);
                    end
                else
                    if isstruct(in2) %parse input structure
                        p.StructExpand = true;
                        p.KeepUnmatched = true;
                    end
                    for field_index = 1:length(field_names)
                        field = field_names{field_index};
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
        end % function parse_input()
        
        
        function out_param = validate_inputs(out_param)
            % let end point of interval not be infinity
            default = out_param.default_values;
            hyperbox = out_param.hyperbox;
            if any(isnan(hyperbox(:))); %check hyperbox for not a number
                out_param.exit=10; out_param = cubMC_g_err(out_param); return;
            end
            [two, out_param.dim]=size(hyperbox); %hyperbox should be 2 x dimension
            if two==0 && isfield(out_param,'hyperbox');
                %if hyperbox specified through out_param structure
                hyperbox=out_param.hyperbox; %then get it from there
                [two, out_param.dim]=size(hyperbox); %and get the dimension
            end
            if two~=2 %if hyperbox is given as row vector for dimension 1, fix that
                if out_param.dim==2; out_param.dim=two; hyperbox=hyperbox';
                else out_param.exit=11; out_param = cubMC_g_err(out_param); return;
                    %else, return an error
                end
            end
            hyperbox=[min(hyperbox,[],1); max(hyperbox,[],1)];
            %ensure left and right endpoints are in order
            if any(hyperbox(1,:)==hyperbox(2,:)); %hyperbox is a point in one direction
                out_param.exit=12; out_param = cubMC_g_err(out_param); return;
            end
            out_param.hyperbox=hyperbox; %copy hyperbox into the out_param structure
            
            if isfield(out_param,'measure'); % the sample measure
                out_param.measure = validatestring(out_param.measure,{'uniform','normal','Gaussian'});
                if strcmpi(out_param.measure,'Gaussian')
                    out_param.measure='normal';
                end
            else
                out_param.measure=default.measure;
            end
            if strcmpi(out_param.measure,'uniform')&&~all(isfinite(hyperbox(:)))
                %cannot integrate on an infinite hyperbox with the uniform distribution
                out_param.exit=13; out_param = cubMC_g_err(out_param); return;
            end
            if strcmpi(out_param.measure,'normal')&&any(isfinite(hyperbox(:)))
                %must integrate on an infinite hyperbox with the normal distribution
                out_param.exit=14; out_param = cubMC_g_err(out_param); return;
            end
            
            if (out_param.abstol < 0)
                %absolute error tolerance should be positive
                warning('GAIL:cubMC_g:abstolneg',...
                    ['The absolute error tolerance should be positive; '...
                    'We will take the absolute value of the absolute error tolerance provided.'])
                out_param.abstol = abs(out_param.abstol);
            end
            if (out_param.reltol < 0 || out_param.reltol > 1)
                % relative error tolerance should be in [0,1]
                warning('GAIL:cubMC_g:reltolneg',...
                    ['Relative error tolerance should be in [0,1]; ' ...
                    'We will use the default value of the error tolerance 1e-1.'])
                out_param.reltol = default.reltol;
            end
            if (out_param.fudge <= 1)
                %standard deviation inflation factor should be a number bigger than 1
                warning('GAIL:cubMC_g:fudgelessthan1',...
                    ['The fudge factor should be bigger than 1; '...
                    'We will use the default value 1.2.'])
                out_param.fudge = default.fudge;
            end
        end % function validate_inputs
        
        
        function out_struct = toStruct(out_param,varargin)
            field_list = ...%union({'f'}, union(out_param.input_field_names, out_param.output_field_names,'stable'),'stable');
                {'f','measure','abstol','reltol','fudge','dim', 'hyperbox'
                };
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
        
         function default = get_default(out_param)
            default = out_param.default_values;
        end
        function out_param = set_input_field_names(out_param, names)
            out_param.input_field_names = names;
        end
        

    end % methods
    

    methods (Access = protected) % seen by subclasses
        
        % customize display order of properties (data) in an instance
        function propgrp = getPropertyGroups(~)
            proplist = {'f', 'measure','abstol','reltol','fudge','dim','hyperbox'};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end

    end % methods (protected)

    
end % classdef

%% other functions