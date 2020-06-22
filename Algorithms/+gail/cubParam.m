classdef cubParam < handle & matlab.mixin.CustomDisplay
   %GAIL.CUBPARAM is a class containing the parameters related to
   %algorithms that find the mean of a random variable
   %   This class contains the number of integrands with the same integral,
   %   etc.
   %
   % Example 1. Construct a cubParam object with default parameters
   % >> cubParamObj = gail.cubParam
   % cubParamObj = ***
   %
   %           f: @(x)sum(x.^2,2)
   %      domain: [2***1 double]
   % measureType: 'uniform'
   %     measure: 'uniform'
   %      absTol: 0.0100
   %      relTol: 0
   %
   %
   % Example 2. Construct a cubParam object with properly ordered inputs
   % >> cubParamObj = gail.cubParam(@(x) sum(x.^3.2),[0 0; 2 2],'box','Lebesgue')
   % cubParamObj = ***
   %
   %              f: @(x)sum(x.^3.2)
   %         domain: [2***2 double]
   %    measureType: 'uniform'
   %        measure: 'Lebesgue'
   %         absTol: 0.0100
   %         relTol: 0
   %
   %
   % Example 3. Using name/value pairs
   % >> cubParamObj = gail.cubParam('domain', [-Inf -Inf; Inf Inf], 'f', @(x) sum(x.^3.2), 'relTol', 0.1, 'measure', 'Gaussian')
   % cubParamObj = ***
   %
   %              f: @(x)sum(x.^3.2)
   %         domain: [2***2 double]
   %    measureType: 'uniform'
   %        measure: 'normal'
   %         absTol: 0.0100
   %         relTol: 0.1000
   %
   %
   % Example 4. Using a structure for input
   % >> inpStruct.f = @(x) sin(sum(x,2));
   % >> inpStruct.domain = [zeros(1,4); ones(1,4)];
   % >> inpStruct.nInit = 2048;
   % >> cubParamObj = gail.cubParam(inpStruct)
   % cubParamObj = ***
   %
   %              f: @(x)sin(sum(x,2))
   %         domain: [2***4 double]
   %    measureType: 'uniform'
   %        measure: 'uniform'
   %         absTol: 0.0100
   %         relTol: 0
   %          nInit: 2048
   %
   %
   % Example 5. Copying a cubParam object and changing some properties
   % >> NewCubParamObj = gail.cubParam(cubParamObj,'measure','Lebesgue')
   % NewCubParamObj = ***
   %
   %              f: @(x)sin(sum(x,2))
   %         domain: [2***4 double]
   %    measureType: 'uniform'
   %        measure: 'Lebesgue'
   %         absTol: 0.0100
   %         relTol: 0
   %          nInit: 2048
   %
   %
   % Author: Fred J. Hickernell

   properties
      err %an errorParam object
      fun %an fParam object
      CM %a cubMeanParam object
      measure %measure against which to integrate
      measureType %type of measure that final integral is made with respect to
      nf %number of f for each integral
      exitflag % exit flag
   end

   properties (Dependent = true)
      ff %function after variable transformation
      nCV %number of control variates
      volume %volume of the domain with respect to the measure
      mmin
      mmax
   end

   properties (Hidden, SetAccess = private)
      def_measure = 'uniform' % default measure
      def_measureType = 'uniform' % default measure
      def_nf = 1 %default number of functions per integral
      allowedMeasures = {'uniform',  ... %over a hyperbox volume of domain is one
         'Lebesgue', ... %like uniform, but integral over domain is the volume of the domain
         'Gaussian', 'normal' ... %these are the same
         }
      allowedmeasureTypes = {'uniform',  ... %over a hyperbox volume of domain is one
         'Gaussian', 'normal' ... %these are the same
         }
   end

   methods
      % Creating a cubParam process
      function obj = cubParam(varargin)
         %this constructor essentially parses inputs
         %the parser will look for the following in order
         %  # a copy of a cubParam object
         %  # a structure
         %  # a function
         %  # a domain
         %  # a domain type
         %  # a measure
         %  # numbers: absTol, relTol, trueMuCV, nMu, nf, inflate
         %  # name-value pairs

         start = 1; %index to begin to parse
         useDefaults = true; %true unless copying an fParam object, then false
         objInp = 0; %where is the an object in the class
         fInp = 0; %where is the input function
         domainInp = 0; %where is the domain
         domainTypeInp = 0; %where is the domain type
         measureInp = 0; %where is the integration measure
         structInp = 0; %where is the structure
         if nargin %there are inputs to parse and assign
            if isa(varargin{start},'gail.cubParam')
               %the first input is a meanYParam object so copy it
               objInp = start;
               start = start + 1;
            end
            if nargin >= start
               if isstruct(varargin{start}) %next input is a structure containing Y
                  structInp = start;
                  start = start + 1;
               end
            end
            if nargin >= start
               if gail.isfcn(varargin{start}) %next input is the function f
                  fInp = start;
                  start = start + 1;
                  if nargin >= start
                     if ismatrix(varargin{start}) && numel(varargin{start}) > 1
                        %next input is the domain
                        domainInp = start;
                        start = start+1;
                        if nargin >= start
                           if ischar(varargin{start}) %next imput is domain type
                              domainTypeInp = start;
                              start = start+ 1;
                              if nargin >= start
                                 if ischar(varargin{start}) %next input is meaure
                                    measureInp = start;
                                    start = start + 1;
                                 end
                              end
                           end
                        end
                     end
                  end
               end
            end
         end

         if objInp
            val = varargin{objInp}; %first input
            obj.err = gail.errorParam(val.err);
            obj.fun = gail.fParam(val.fun);
            obj.CM = gail.cubMeanParam(val.CM);
            obj.measure = val.measure; %copy integration measure
            obj.measureType = val.measureType; %copy integration measure
            obj.nf = val.nf; %copy number of functions for each integral
            obj.exitflag=val.exitflag;
            useDefaults = false;
         end

         %Parse fParam properties
         whichfParse = [fInp domainInp domainTypeInp structInp start:nargin];
         whichfParse = whichfParse(whichfParse > 0);
         if objInp
            obj.fun = gail.fParam(obj.fun,varargin{whichfParse});
         else
            obj.fun = gail.fParam(varargin{whichfParse});
         end

         %Parse errorParam properties
         whichErrParse = [structInp start:nargin];
         whichErrParse = whichErrParse(whichErrParse > 0);
         if objInp
            obj.err = gail.errorParam(obj.err,varargin{whichErrParse});
         else
            obj.err = gail.errorParam(varargin{whichErrParse});
         end

         %Parse cubMeanParam properties
         whichCMParse = [structInp start:nargin];
         whichCMParse = whichCMParse(whichCMParse > 0);
         if objInp
            obj.CM = gail.cubMeanParam(obj.CM,varargin{whichCMParse});
         else
            obj.CM = gail.cubMeanParam(varargin{whichCMParse});
         end

         %Now begin to parse remaining inputs
         p = inputParser; %construct an inputParser object
         p.KeepUnmatched = true; %ignore those that do not match
         p.PartialMatching = false; %don'try a partial match
         p.StructExpand = true; %expand structures
         done = false; %not finished parsing
         if nargin >= start
            if ischar(varargin{start})
               %there may be input string/value pairs or a structure
               MATLABVERSION = gail.matlab_version;
               if MATLABVERSION >= 8.3
                  f_addParamVal = @addParameter;
               else
                  f_addParamVal = @addParamValue;
               end
               parseRange = start:nargin;
               done = true;
            end
         end

         if ~done %then nothingleft or just numbers
            f_addParamVal = @addOptional;
            parseRange = []; %nothing to parse here if just numbers
         end

         f_addParamVal(p,'measure',obj.def_measure);
         f_addParamVal(p,'measureType',obj.def_measureType);
         f_addParamVal(p,'nf',obj.def_nf);

         if structInp
            parse(p,varargin{parseRange},varargin{structInp})
            %parse inputs with a structure
         else
            parse(p,varargin{parseRange}) %parse inputs w/o structure
         end
         struct_val = p.Results; %store parse inputs as a structure
         if ~useDefaults %remove defaults if copying cubParam object
            struct_val = rmfield(struct_val,p.UsingDefaults);
         end

         %Assign values of structure to corresponding class properties
         if measureInp
            obj.measure = varargin{measureInp}; %assign measure
         elseif isfield(struct_val,'measure')
            obj.measure = struct_val.measure;
         end
         if isfield(struct_val,'measureType')
            obj.measureType = struct_val.measureType;
         end
         if isfield(struct_val,'nf')
            obj.nf = struct_val.nf;
         end

      end %of constructor

      function set.measure(obj,val)
         validateattributes(val, {'char'}, {})
         obj.measure = checkMeasure(obj,val);
      end

      function set.measureType(obj,val)
         validateattributes(val, {'char'}, {})
         obj.measureType = checkmeasureType(obj,val);
      end

      function set.nf(obj,val)
         validateattributes(val, {'numeric'}, {'positive','integer'})
         obj.nf = checkNf(obj,val);
      end

      function val = get.volume(obj) %volume of the domain
         if any(strcmp(obj.measure,{'uniform', 'normal'}))
            if strcmpi(obj.fun.domainType, 'ball')
               val = ((2.0*pi^(obj.d/2.0))/(obj.d*gamma(obj.d/2.0)))*obj.radius^obj.d; %volume of a d-dimentional ball
            elseif strcmpi(obj.fun.domainType, 'sphere')
               val = ((2.0*pi^(obj.d/2.0))/(gamma(obj.d/2.0)))*obj.radius^(obj.d - 1); %volume of a d-dimentional sphere
            end
         elseif strcmp(obj.measure, {'Lebesgue'})
            val = prod(diff(obj.fun.domain,1),2);
         end
      end

      % Mmin and mmax is based on nInit and nMax
      function val=get.mmin(obj)
         val = ceil(log2(obj.CM.nInit));
      end

      function val=get.mmax(obj)
         val=floor(log2(obj.CM.nMax));
      end

      % Transformation section
      function val = get.ff(obj)

         % Domain type = box
         if strcmp(obj.fun.domainType,'box')
            if strcmp(obj.measure,'uniform')
               Cnorm = prod(obj.fun.domain(2,:)-obj.fun.domain(1,:));
               val =@(t) Cnorm*obj.fun.f(bsxfun(@plus,obj.fun.domain(1,:), ...
                  bsxfun(@times,(obj.fun.domain(2,:)-obj.fun.domain(1,:)),t))); % a + (b-a)x = u

            elseif strcmp(obj.measure,'Lebesgue')
               Cnorm = prod(obj.fun.domain(2,:)-obj.fun.domain(1,:));
               val =@(t) obj.volume*Cnorm*obj.fun.f(bsxfun(@plus,obj.fun.domain(1,:), ...
                  bsxfun(@times,(obj.fun.domain(2,:)-obj.fun.domain(1,:)),t))); % a + (b-a)x = u

            elseif strcmp(obj.measure, 'normal')
               val = @(t) obj.fun.f(gail.stdnorminv(t));
            end
         end

         % Domain type ball/sphere - automatically assume the original
         % measure is uniform. If it is lesbegue, then we will just
         % multiple by volume at next step
         if strcmp(obj.measureType, 'uniform')
            if strcmp(obj.fun.domainType, 'ball') || strcmp(obj.fun.domainType, 'ball-from-cube')
               val=@(t) obj.fun.f(gail.domain_balls_spheres.ball_psi_1(t, obj.d, obj.radius, obj.fun.domain));

            elseif strcmp(obj.fun.domainType, 'ball-from-normal')
               val=@(t) obj.fun.f(gail.domain_balls_spheres.ball_psi_2(gail.stdnorminv(t), obj.d, obj.radius, obj.fun.domain));

            elseif strcmp(obj.fun.domainType,'sphere') || strcmp(obj.fun.domainType, 'sphere-from-cube')
               val=@(t) obj.fun.f(gail.domain_balls_spheres.sphere_psi_1(t, obj.d, obj.radius, obj.fun.domain));

            elseif strcmp(obj.fun.domainType, 'sphere-from-normal')
               val=@(t) obj.fun.f(gail.domain_balls_spheres.sphere_psi_2(gail.stdnorminv(t), obj.d, obj.radius, obj.fun.domain));
            end

         elseif strcmp(obj.measureType, 'normal')
            if strcmp(obj.fun.domainType, 'ball') || strcmp(obj.fun.domainType, 'ball-from-normal')
               val=@(t) obj.fun.f(gail.domain_balls_spheres.ball_psi_2(t, obj.d, obj.radius, obj.fun.domain));

            elseif strcmp(obj.fun.domainType,'ball-from-cube')
               val=@(t) obj.fun.f(gail.domain_balls_spheres.ball_psi_2(gail.stdnorm(t), obj.d, obj.radius, obj.fun.domain));

            elseif strcmp(obj.fun.domainType,'sphere') || strcmp(obj.fun.domainType, 'sphere-from-normal')
               val=@(t) obj.fun.f(gail.domain_balls_spheres.sphere_psi_2(t, obj.d, obj.radius, obj.fun.domain));

            elseif strcmp(obj.fun.domainType, 'sphere-from-cube')
               val=@(t) obj.fun.f(gail.domain_balls_spheres.sphere_psi_2(gail.stdnorm(t), obj.d, obj.radius, obj.fun.domain));
            end
         end

         % If original measure if lesbeque, just multiple by volume
         if strcmp(obj.measure, 'Lesbeque')
            val = obj.volume*val;

         elseif strcmp(obj.measure, 'normal')
            val = @(t) obj.fun.f(gail.stdnorminv(t));
         end

      end
   end

   methods (Access = protected)
      function outval = checkMeasure(obj,inval)
         assert(any(strcmp(inval,obj.allowedMeasures)))
         if strcmp(inval,'Gaussian') %same as normal
            outval = 'normal';
         else
            outval = inval;
         end

         if strcmp(outval,'normal') %domain must be R^d
            assert(all(obj.fun.domain(1,:) == -Inf) && ...
               all(obj.fun.domain(2,:) == Inf))
         elseif any(strcmp(outval,{'uniform','Lebesgue'})) %domain must be finite
            assert(all(all(isfinite(obj.fun.domain))))
         end
      end

      function outval = checkmeasureType(obj,inval)
         assert(any(strcmp(inval,obj.allowedmeasureTypes)))
         if strcmp(inval,'Gaussian') %same as normal
            outval = 'normal';
         else
            outval = inval;
         end
      end

      function outval = checkNf(obj,inval)
         assert(obj.fun.nfOut - sum(inval) == obj.CM.nCV)
         outval = inval;
      end

      function propgrp = getPropertyGroups(obj)
         if ~isscalar(obj)
            propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
         else
            propList = getPropertyList(obj);
            propgrp = matlab.mixin.util.PropertyGroup(propList);
         end
      end

      function propList = getPropertyList(obj)
         propList = struct('f', obj.fun.f, ...
            'domain', obj.fun.domain);
         if ~strcmp(obj.fun.domainType,obj.fun.def_domainType)
            propList.domainType = obj.fun.domainType;
         end

         propList.measureType=obj.measureType;
         propList.measure = obj.measure;
         propList.absTol = obj.err.absTol;
         propList.relTol = obj.err.relTol;

         if obj.CM.nInit ~= obj.CM.def_nInit
            propList.nInit = obj.CM.nInit;
         end
         if obj.CM.nMax ~= obj.CM.def_nMax
            propList.nMax = obj.CM.nMax;
         end
         if obj.CM.nMu ~= obj.CM.def_nMu
            propList.nMu = obj.CM.nMu;
         end
         if obj.nf ~= obj.def_nf
            propList.nf= obj.nf;
         end
         if numel(obj.CM.trueMuCV)
            propList.trueMuCV = obj.CM.trueMuCV;
         end
      end

   end

end
