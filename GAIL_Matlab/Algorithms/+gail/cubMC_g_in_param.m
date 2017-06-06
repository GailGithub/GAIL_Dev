% cubMC_g_in_param: cubMC_g's input parameter object
%
% Examples
%
% >> in_param = gail.cubMC_g_in_param()
%     Warning: Function f must be a function handle. Now GAIL is using f(x)=exp(-100*(x-0.5)^2). 
%     ***
% 
%     in_param = 
% 
%       cubMC_g_in_param with properties:
% 
%          measure: 'uniform'
%           abstol: 0.0100
%           reltol: 0.1000
%            alpha: 0.0100
%            fudge: 1.2000
%             nSig: 10000
%               n1: 10000
%          tbudget: 100
%          nbudget: 1.0000e+09
%             flag: 2
%              dim: 1
%         hyperbox: [2x1 double]
%
%
%
% >> f = @(x) x.^2; in_param = gail.cubMC_g_in_param(f)
%    in_param =***
%          measure: 'uniform'
%           abstol: 0.0100
%           reltol: 0.1000
%            alpha: 0.0100
%            fudge: 1.2000
%             nSig: 10000
%               n1: 10000
%          tbudget: 100
%          nbudget: 1.0000e+09
%             flag: 2
%              dim: 1
%         hyperbox: [2x1 double]
%
%
%  To get a struct:
%  >> in_param = gail.cubMC_g_in_param(@(x) x.^2);  out_param = in_param.toStruct()
%  out_param =***
%          measure: 'uniform'
%           abstol: 0.0100
%           reltol: 0.1000
%            alpha: 0.0100
%            fudge: 1.2000
%             nSig: 10000
%               n1: 10000
%          tbudget: 100
%          nbudget: 1.0000e+09
%             flag: 2
%              dim: 1
%         hyperbox: [2x1 double]
%
%
% To get a structure with selected fields (and ignore properties that do not exist):
% >> in_param = gail.cubMC_g_in_param(@(x) x.^2);  out_param = in_param.toStruct({'f','measure','hyperbox','nonexistent'})
%  out_param =***
%            f: @(x)x.^2
%      measure: 'uniform'
%     hyperbox: [2x1 double]
%
%
% >> f=@(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [0 0; 1 1];
% >> in_param = gail.cubMC_g_in_param(f,hyperbox,'measure','uniform','abstol',1e-3,'reltol',0)
%    in_param =***
%          measure: 'uniform'
%           abstol: 1.0000e-03
%           reltol: 0
%            alpha: 0.0100
%            fudge: 1.2000
%             nSig: 10000
%               n1: 10000
%          tbudget: 100
%          nbudget: 1.0000e+09
%             flag: 2
%              dim: 2
%         hyperbox: [2x2 double]
%    
%
% >> f=@(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-inf -inf;inf inf];
% >> in_param = gail.cubMC_g_in_param(f,hyperbox,'normal',0,1e-2)
%    in_param =***
%          measure: 'normal'
%           abstol: 0
%           reltol: 0.0100
%            alpha: 0.0100
%            fudge: 1.2000
%             nSig: 10000
%               n1: 10000
%          tbudget: 100
%          nbudget: 1.0000e+09
%             flag: 2
%              dim: 2
%         hyperbox: [2x2 double]
%
%    >> in_param.hyperbox
% 
%     ans =
% 
%       -Inf  -Inf
%        Inf   Inf
%
% >> d=3;f=@(x) 2^d*prod(x,2)+0.555; hyperbox =[zeros(1,d);ones(1,d)];
% >> in_param.abstol = 1e-3; in_param.reltol=1e-3;
% >> in_param = gail.cubMC_g_in_param(f,hyperbox,in_param)
%          measure: 'uniform'
%           abstol: 0.0100
%           reltol: 0.1000
%            alpha: 0.0100
%            fudge: 1.2000
%             nSig: 10000
%               n1: 10000
%          tbudget: 100
%          nbudget: 1.0000e+09
%             flag: 2
%              dim: 3
%         hyperbox: [2x3 double]
%
%    >> in_param.hyperbox
% 
%     ans =
% 
%         0     0     0
%         1     1     1
%
classdef cubMC_g_in_param < gail.gailMD_in_param
    %% data
    properties % public
        alpha
        nSig
        n1
        tbudget
        nbudget
        flag
    end % properties
     
    %% methods
    methods % public
        % constructor
        function out_param = cubMC_g_in_param(varargin)
            % parse the input to a gail function
            in = cell(0);
            if nargin >= 1
                in = varargin{1};
            end
            out_param = out_param@gail.gailMD_in_param(in);

            % out_param.get_input_field_names();
            %% Default parameter values
            default = out_param.get_default();
            default.alpha = 0.01;% default uncertainty
            default.nSig = 1e4; % default nSig initial sample size to estimate sigma
            default.n1 = 1e4; % default n1 initial sample size to estimate Q
            default.tbudget = 100;% default time budget in seconds
            default.nbudget = 1e9; % default sample budget
            default.flag = 0; % default value of parameter checking status
        
            %% parse inputs
            out_param = out_param.set_input_field_names(...
              {'measure','abstol','reltol','alpha','fudge','nSig','n1','tbudget','nbudget','flag','dim'}...
            );
            out_param = out_param.parse_inputs(default, varargin{:});
            
            %% validate inputs

            out_param = out_param.validate_inputs();
            
        end % constructor
        
        
        function out_param = toStruct(out_param,varargin)
            l = {'measure','abstol','reltol','alpha','fudge','nSig','n1','tbudget','nbudget','flag','dim','hyperbox'};
            if ~isempty(varargin)   
                l = varargin{1};
            end
            out_param = toStruct@gail.gailMD_in_param(out_param, l);
        end
        
        function out_param = validate_inputs(out_param)
            
            out_param = validate_inputs@gail.gailMD_in_param(out_param);
            
            if (~isempty(out_param.alpha))
                if out_param.flag == 0
                    if (out_param.alpha <= 0 || out_param.alpha >= 1)
                        %uncertainty should be 1 (0,1)
                        warning('GAIL:cubMC_g:alphanot01',...
                            ['The uncertainty should be in (0,1); '...
                            'We will use the default value 1e-2.'])
                        out_param.alpha = default.alpha;
                    end
                    if (~gail.isposge30(out_param.nSig))
                        %the sample to estimate sigma should be a positive integer
                        warning('GAIL:cubMC_g:nsignotposint',...
                            ['The number nSig should a positive integer greater than 30; '...
                            'We will take the default value 1e4.'])
                        out_param.nSig = default.nSig;
                    end
                    if (~gail.isposge30(out_param.n1))
                        %initial sample size to estimate Q should be a positive integer
                        warning('GAIL:cubMC_g:n1notposint',...
                            ['The number n1 should a positive integer greater than 30; '...
                            'We will use the default value 1e4.'])
                        out_param.n1 = default.n1;
                    end
                    if (out_param.tbudget < 0)
                        %the time budget in seconds should be positive
                        warning('GAIL:cubMC_g:tbudgetneg',...
                            ['Time budget should be positive; '...
                            'We will use the absolute value of the time budget'])
                        out_param.tbudget = abs(out_param.tbudget);
                    end
                    if (out_param.tbudget == 0)
                        %the time budget in seconds should be positive
                        warning('GAIL:cubMC_g:tbudget0',...
                            ['Time budget should be positive rather than 0; '...
                            'We will use the default value of the time budget 100 seconds.'])
                        out_param.tbudget = default.tbudget;
                    end
                    if (~gail.isposge30(out_param.nbudget))
                        %the sample budget should be a positive integer
                        warning('GAIL:cubMC_g:nbudgetnotposint',...
                            ['The number of sample budget should be a large positive integer;'...
                            'We will take the default value of 1e9.'])
                        out_param.nbudget = default.nbudget;
                    end
                end
                out_param.flag = 2;
            end
        end

    end % methods
    
    methods (Access = protected) % seen by subclasses
        
        % customize display order of properties (data) in an instance
        function propgrp = getPropertyGroups(~)
            proplist =  {'measure','abstol','reltol','alpha','fudge','nSig','n1','tbudget','nbudget','flag','dim','hyperbox'};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end % methods (protected)
    
    
end % classdef

%% other functions