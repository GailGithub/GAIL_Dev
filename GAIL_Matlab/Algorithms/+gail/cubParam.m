classdef cubParam < gail.fParam
   %GAIL.CUBPARAM is a class containing the parameters related to
   %algorithms that find the mean of a random variable
   %   This class contains the number of integrands with the same integral,
   %   etc.
   %
   % Example 1. Construct a cubParam object with default parameters
   % >> cubParamObj = gail.cubParam
   % cubParamObj =
   %   cubParam with properties:
   %
   %              f: @(x)sum(x.^2,2)
   %         domain: [2×1 double]
   %        measure: 'uniform'
   %         absTol: 0.010000000000000
   %         relTol: 0
   %
   %
   % Example 2. Construct a cubParam object with properly ordered inputs
   % >> cubParamObj = gail.cubParam(@(x) sum(x.^3.2),[0 0; 2 2],'box','Lebesgue')
   % cubParamObj =
   %   cubParam with properties:
   %
   %              f: @(x)sum(x.^3.2)
   %         domain: [2×2 double]
   %        measure: 'Lebesgue'
   %         absTol: 0.010000000000000
   %         relTol: 0
   %
   %
   % Example 3. Using name/value pairs
   % >> cubParamObj = gail.cubParam('domain', [-Inf -Inf; Inf Inf], 'f', @(x) sum(x.^3.2), 'relTol', 0.1, 'measure', 'Gaussian')
   % cubParamObj =
   %   cubParam with properties:
   %
   %              f: @(x)sum(x.^3.2)
   %         domain: [2×2 double]
   %        measure: 'normal'
   %         absTol: 0.010000000000000
   %         relTol: 0.100000000000000
   %
   %
   % Example 4. Using a structure for input
   % >> inpStruct.f = @(x) sin(sum(x,2));
   % >> inpStruct.domain = [zeros(1,4); ones(1,4)];
   % >> inpStruct.nInit = 1000;
   % >> cubParamObj = gail.cubParam(inpStruct)
   % cubParamObj =
   %   cubParam with properties:
   %
   %              f: @(x)sin(sum(x,2))
   %         domain: [2×4 double]
   %        measure: 'uniform'
   %         absTol: 0.010000000000000
   %         relTol: 0
   %
   %
   % Example 5. Copying a cubParam object and changing some properties
   % >> NewCubParamObj = gail.cubParam(cubParamObj,'measure','Lebesgue')
   % NewCubParamObj =
   %   cubParam with properties:
   %
   %              f: @(x)sin(sum(x,2))
   %         domain: [2×4 double]
   %        measure: 'Lebesgue'
   %         absTol: 0.010000000000000
   %         relTol: 0
   %          nInit: 1000
   %
   %
   
   
   % Author: Fred J. Hickernell
   
   properties
      measure %measure against which to integrate
      trueMuCV %true integral for control variates
      nMu %number of integrals for solution function
      nf %number of f for each integral
      inflate %inflation factor for bounding the error
      ffMeasure %measure for the transformed integrand
      
      % Willy added
      radius
      transf

      betaUpdate
      
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
      def_inflate = @(m) 5 * 2^-m %default inflation factor
      def_nMu = 1 %default number of integrals
      def_nf = 1 %default number of Y per integral
      def_trueMuCV = [] %default true integrals for control variates
      def_ffMeasure ='uniform'
      def_radius = 1
      def_transf = 1
      
      % Willy added
      def_betaUpdate=0
      
      allowedMeasures = {'uniform', 'uniform ball', ... %over a hyperbox volume of domain is one
         'Lebesgue', ... %like uniform, but integral over domain is the volume of the domain
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
         

         %Parse errorParam properties
         whichParse = [objInp fInp domainInp domainTypeInp structInp start:nargin];
         whichParse = whichParse(whichParse > 0);
         obj@gail.fParam(varargin{whichParse});
         
         if objInp
            val = varargin{objInp}; %first input
            obj.measure = val.measure; %copy integration measure
            obj.inflate = val.inflate; %copy inflation factor
            obj.nMu = val.nMu; %copy the number of integrals
            obj.nf = val.nf; %copy number of functions for each integral
            obj.trueMuCV = val.trueMuCV; %copy true means of control variates
            obj.ffMeasure= val.ffMeasure;
            
            % Willy added
            obj.betaUpdate=val.betaUpdate;
            useDefaults = false;
            
         end
         
         %Now begin to parse inputs
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
            parseRange = (start + 2):nargin; %to account for the two tolerances already parsed
         end
         
         f_addParamVal(p,'measure',obj.def_measure);
         f_addParamVal(p,'trueMuCV',obj.def_trueMuCV);
         f_addParamVal(p,'nMu',obj.def_nMu);
         f_addParamVal(p,'nf',obj.def_nf);
         f_addParamVal(p,'inflate',obj.def_inflate);
         f_addParamVal(p,'ffMeasure',obj.def_ffMeasure);

         % Willy added
         f_addParamVal(p,'radius',obj.def_radius);
         f_addParamVal(p,'transf',obj.def_transf);
         f_addParamVal(p,'betaUpdate', obj.def_betaUpdate);
         

         
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
         if isfield(struct_val,'inflate')
            obj.inflate = struct_val.inflate;
         end
         if isfield(struct_val,'nMu')
            obj.nMu = struct_val.nMu;
         end
         if isfield(struct_val,'nf')
            obj.nf = struct_val.nf;
         end
         
         if isfield(struct_val,'trueMuCV')
            obj.trueMuCV = struct_val.trueMuCV;
         end
        
         if isfield(struct_val,'ffMeasure')
            obj.ffMeasure = struct_val.ffMeasure;
         end
         
         % Willy added
         if isfield(struct_val,'radius')
            obj.radius = struct_val.radius;
         end
         
         if isfield(struct_val,'transf')
            obj.transf = struct_val.transf;
         end
         
         if isfield(struct_val,'betaUpdate')
            obj.betaUpdate = struct_val.betaUpdate;
         end
         
      end %of constructor
      
      function set.measure(obj,val)
         validateattributes(val, {'char'}, {})
         obj.measure = checkMeasure(obj,val);
      end
      
      function set.ffMeasure(obj,val)
         validateattributes(val, {'char'}, {})
         obj.ffMeasure = val;
      end
      
      function set.inflate(obj,val)
         validateattributes(val, {'function_handle'}, {})
         obj.inflate = val;
      end
      
      function set.nMu(obj,val)
         validateattributes(val, {'numeric'}, {'integer', 'positive'})
         obj.nMu = val;
      end
      
      function set.nf(obj,val)
         validateattributes(val, {'numeric'}, {'positive','integer'})
         obj.nf = val;
      end
      
      
      function set.trueMuCV(obj,val)
         validateattributes(val, {'numeric'}, {})
         obj.trueMuCV = setTrueMuCVDim(obj,val);
      end
      
      function set.betaUpdate(obj,val)
         validateattributes(val, {'numeric'}, {})
         obj.betaUpdate = val;
      end
      
      function val = get.nCV(obj)
         val = obj.nfOut - sum(obj.nf);
      end
      
      function val=get.mmin(obj)
         val=ceil(log2(obj.nInit));
      end 
      
      function val=get.mmax(obj)
         val=floor(log2(obj.nMax));
      end 
         
      function val = get.volume(obj) %volume of the domain
         
         %% Working with Fred 7/10/2017
%          if any(strcmp(obj.measure,{'uniform', 'normal'}))
%             if strcmpi(obj.domainType, 'ball')
%                val = ((2.0*pi^(obj.d/2.0))/(obj.d*gamma(obj.d/2.0)))*obj.radius^obj.d; %volume of a d-dimentional ball
%             elseif strcmpi(obj.domainType, 'sphere')
%                val = ((2.0*pi^(obj.d/2.0))/(gamma(obj.d/2.0)))*obj.radius^(obj.d - 1); %volume of a d-dimentional sphere
%             end
%             
%          elseif strcmp(obj.measure, {'Lebesgue'})
%             val = prod(diff(obj.domain,1),2);
%          end
         
         %          ORIGINAL WORKING (7/8 Examples)
%                   if any(strcmp(obj.measure,{'uniform', 'normal'}))
%                      val = 1;
%                   elseif strcmp(obj.measure, {'Lebesgue'})
%                      val = prod(diff(obj.domain,1),2);
%                   end
                
                  if strcmpi(obj.measure,'uniform ball')% using the formula of the obj.volume of a ball
                    val = ((2.0*pi^(obj.d/2.0))/(obj.d*gamma(obj.d/2.0)))*obj.radius^obj.d; %obj.volume of a d-dimentional ball
                  else % using the formula of the obj.volume of a sphere
                    val = ((2.0*pi^(obj.d/2.0))/(gamma(obj.d/2.0)))*obj.radius^(obj.d - 1); %obj.volume of a d-dimentional sphere
                  end
         
      end
      
      function val = get.ff(obj)
         
         %% Working with Fred 7/10/2017
         %          if strcmp(obj.domainType,'box')
         %             if strcmp(obj.measure,'uniform')
         %                val = @(t) obj.f(bsxfun(@plus, obj.domain(1,:), ...
         %                   bsxfun(@times, diff(obj.domain,1), t)));
         %             elseif strcmp(obj.measure,'Lebesgue')
         %                val = @(t) obj.volume*obj.f(bsxfun(@plus, obj.domain(1,:), ...
         %                   bsxfun(@times, diff(obj.domain,1), t)));
         %             elseif strcmp(obj.measure, 'normal')
         %                val = @(t) obj.f(gail.stdnorminv(t));
         %             end
         %          end
         
%          if strcmp(obj.measure, 'uniform')
%             % Box section
%             if strcmp(obj.domainType, 'cube')
%                Cnorm = prod(obj.domain(2,:)-obj.domain(1,:));
%                val =@(t) Cnorm*obj.f(bsxfun(@plus,obj.domain(1,:), ...
%                   bsxfun(@times,(obj.domain(2,:)-obj.domain(1,:)),t))); % a + (b-a)x = u
% 
%                
%             end
%             
%             % ffMeasure uniform
%             if strcmp(obj.ffMeasure, 'uniform')
%                % Ball section
%                if strcmp(obj.domainType, 'ball')
%                   val=@(t) obj.f(gail.domain_balls_spheres.ball_psi_1(t, obj.d, obj.radius, obj.domain));
%                elseif strcmp(obj.domainType, 'ball-from-cube')
%                   val=@(t) obj.f(gail.domain_balls_spheres.ball_psi_1(t, obj.d, obj.radius, obj.domain));
%                elseif strcmp(obj.domainType, 'ball-from-normal')
%                   val=@(t) obj.f(gail.domain_balls_spheres.ball_psi_2(gail.stdnorminv(t), obj.d, obj.radius, obj.domain));
%                   
%                   % Sphere section
%                elseif strcmp(obj.domainType, 'sphere')
%                   val=@(t) obj.f(gail.domain_balls_spheres.sphere_psi_1(t, obj.d, obj.radius, obj.domain));
%                elseif strcmp(obj.domainType, 'sphere-from-cube')
%                   val=@(t) obj.f(gail.domain_balls_spheres.sphere_psi_1(t, obj.d, obj.radius, obj.domain));
%                elseif strcmp(obj.domainType, 'sphere-from-normal')
%                   val=@(t) obj.f(gail.domain_balls_spheres.sphere_psi_2(gail.stdnorminv(t), obj.d, obj.radius, obj.domain));
%                end
%                
%                % ff measure normal
%             elseif strcmp(obj.ffMeasure, 'normal')
%                % Ball section
%                if strcmp(obj.domainType, 'ball')
%                   val=@(t) obj.f(gail.domain_balls_spheres.ball_psi_2(t, obj.d, obj.radius, obj.domain));
%                elseif strcmp(obj.domainType, 'ball-from-cube')
%                   val=@(t) obj.f(gail.domain_balls_spheres.ball_psi_2(gail.stdnorm(t), obj.d, obj.radius, obj.domain));
%                elseif strmcmp(obj.domainType, 'ball-from-normal')
%                   val=@(t) obj.f(gail.domain_balls_spheres.ball_psi_2(t, obj.d, obj.radius, obj.domain));
%                   
%                   % Sphere section
%                elseif strcmp(obj.domainType, 'sphere')
%                   val=@(t) obj.f(gail.domain_balls_spheres.sphere_psi_2(t, obj.d, obj.radius, obj.domain));
%                elseif strcmp(obj.domainType, 'sphere-from-cube')
%                   val=@(t) obj.f(gail.domain_balls_spheres.sphere_psi_2(gail.stdnorm(t), obj.d, obj.radius, obj.domain));
%                elseif strcmp(obj.domainType, 'sphere-from-normal')
%                   val=@(t) obj.f(gail.domain_balls_spheres.sphere_psi_2(t, obj.d, obj.radius, obj.domain));
%                end
%             end
%             
%          elseif strcmp(obj.measure, 'Lesbesque')
%             % ffMeasure uniform
%             if strcmp(obj.ffMeasure, 'uniform')
%                % Ball section
%                if strcmp(obj.domainType, 'ball')
%                   val=@(t) obj.volume*obj.f(gail.domain_balls_spheres.ball_psi_1(t, obj.d, obj.radius, obj.domain));
%                elseif strcmp(obj.domainType, 'ball-from-cube')
%                   val=@(t) obj.volume*obj.f(gail.domain_balls_spheres.ball_psi_1(t, obj.d, obj.radius, obj.domain));
%                elseif strmcmp(obj.domainType, 'ball-from-normal')
%                   val=@(t) obj.volume*obj.f(gail.domain_balls_spheres.ball_psi_2(gail.stdnorminv(t), obj.d, obj.radius, obj.domain));
%                   
%                   % Sphere section
%                elseif strcmp(obj.domainType, 'sphere')
%                   val=@(t) obj.volume*obj.f(gail.domain_balls_spheres.sphere_psi_1(t, obj.d, obj.radius, obj.domain));
%                elseif strcmp(obj.domainType, 'sphere-from-cube')
%                   val=@(t) obj.volume*obj.f(gail.domain_balls_spheres.sphere_psi_1(t, obj.d, obj.radius, obj.domain));
%                elseif strcmp(obj.domainType, 'sphere-from-normal')
%                   val=@(t) obj.volume*obj.f(gail.domain_balls_spheres.sphere_psi_2(gail.stdnorminv(t), obj.d, obj.radius, obj.domain));
%                end
%                
%                % ff measure normal
%             elseif strcmp(obj.ffMeasure, 'normal')
%                % Ball section
%                if strcmp(obj.domainType, 'ball')
%                   val=@(t) obj.volume*obj.f(gail.domain_balls_spheres.ball_psi_2(t, obj.d, obj.radius, obj.domain));
%                elseif strcmp(obj.domainType, 'ball-from-cube')
%                   val=@(t) obj.volume*obj.f(gail.domain_balls_spheres.ball_psi_2(gail.stdnorm(t), obj.d, obj.radius, obj.domain));
%                elseif strmcmp(obj.domainType, 'ball-from-normal')
%                   val=@(t) obj.volume*obj.f(gail.domain_balls_spheres.ball_psi_2(t, obj.d, obj.radius, obj.domain));
%                   
%                   % Sphere section
%                elseif strcmp(obj.domainType, 'sphere')
%                   val=@(t) obj.volume*obj.f(gail.domain_balls_spheres.sphere_psi_2(t, obj.d, obj.radius, obj.domain));
%                elseif strcmp(obj.domainType, 'sphere-from-cube')
%                   val=@(t) obj.volume*obj.f(gail.domain_balls_spheres.sphere_psi_2(gail.stdnorm(t), obj.d, obj.radius, obj.domain));
%                elseif strcmp(obj.domainType, 'sphere-from-normal')
%                   val=@(t) obj.volume*obj.f(gail.domain_balls_spheres.sphere_psi_2(t, obj.d, obj.radius, obj.domain));
%                end
%             end
%             
%          elseif strcmp(obj.measure, 'normal')
%             val = @(t) obj.f(gail.stdnorminv(t));
%          end
%          

      %         ORIGINAL WORKING (7/8 Examples)
      if strcmpi(obj.measure,'uniform ball') || strcmpi(obj.measure,'uniform sphere') %using uniformly distributed samples on a ball or sphere
         if strcmp(obj.measure,'uniform sphere') && obj.transf == 1 %box-to-sphere transformation
            obj.d = obj.d + 1; % changing obj.d to the dimension of the sphere
            obj.shift = [obj.shift rand];
         end
         
%          if strcmpi(obj.measure,'uniform ball')% using the formula of the obj.volume of a ball
%             obj.volume = ((2.0*pi^(obj.d/2.0))/(obj.d*gamma(obj.d/2.0)))*obj.radius^obj.d; %obj.volume of a d-dimentional ball
%          else % using the formula of the obj.volume of a sphere
%             obj.volume = ((2.0*pi^(obj.d/2.0))/(gamma(obj.d/2.0)))*obj.radius^(obj.d - 1); %obj.volume of a d-dimentional sphere
%          end


         if obj.transf == 1 % box-to-ball or box-to-sphere transformation should be used
            if obj.d == 1 % It is not necessary to multiply the function f by the obj.volume, since no transformation is being made
               obj.domain = [obj.domain - obj.radius; obj.domain + obj.radius];% for one dimension, the ball is actually an interval
               obj.measure = 'uniform';% then a uniform distribution on a box can be used
            else
               if strcmpi(obj.measure,'uniform ball') % box-to-ball transformation
                  val = @(t) obj.f(gail.domain_balls_spheres.ball_psi_1(t, obj.d, obj.radius, obj.domain))*obj.volume;% the psi function is the transformation
               else %  % box-to-sphere transformation
                  val = @(t) obj.f(gail.domain_balls_spheres.sphere_psi_1(t, obj.d, obj.radius, obj.domain))*obj.volume;% the psi function is the transformation
                  obj.d = obj.d - 1;% the box-to-sphere transformation takes points from a (d-1)-dimensional box to a d-dimensional sphere
                  obj.shift = obj.shift(1:end-1);
               end
               obj.domain = [zeros(1, obj.d); ones(1, obj.d)];% the obj.domain must be the domain of the transformation, which is a unit box
               obj.measure = 'uniform';% then a uniform distribution on a box can be used
            end
            
         else % normal-to-ball or normal-to-sphere transformation should be used
            if strcmpi(obj.measure,'uniform ball') % normal-to-ball transformation
               val = @(t) obj.f(gail.domain_balls_spheres.ball_psi_2(t, obj.d, obj.radius, obj.domain))*obj.volume;% the psi function is the transformation
            else % normal-to-sphere transformation
               val = @(t) obj.f(gail.domain_balls_spheres.sphere_psi_2(t, obj.d, obj.radius, obj.domain))*obj.volume;% the psi function is the transformation
            end
            obj.domain = bsxfun(@plus, zeros(2, obj.d), [-inf; inf]);% the obj.domain must be the domain of the transformation, which is a this unit box
            obj.measure = 'normal';% then a normal distribution can be used
         end
      end
            
      if strcmp(obj.measure,'normal')
         val=@(x) obj.f(gail.stdnorminv(x));
      elseif strcmp(obj.measure,'uniform')
         Cnorm = prod(obj.domain(2,:)-obj.domain(1,:));
         val=@(x) Cnorm*obj.f(bsxfun(@plus,obj.domain(1,:),bsxfun(@times,(obj.domain(2,:)-obj.domain(1,:)),x))); % a + (b-a)x = u
      
      
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
            assert(all(obj.domain(1,:) == -Inf) && ...
               all(obj.domain(2,:) == Inf))
         end
         if any(strcmp(outval,{'uniform','Lebesgue'})) %domain must be finite
            assert(all(all(isfinite(obj.domain))))
         end
      end
      
      function outval = setTrueMuCVDim(obj,inval)
         assert(numel(inval) == obj.nCV)
         outval = inval(:)';
      end
      
      function propList = getPropertyList(obj)
         propList = struct('f', obj.f, ...
            'domain', obj.domain);
         if ~strcmp(obj.domainType,obj.def_domainType)
            propList.domainType = obj.domainType;
         end
         
         propList.ffMeasure=obj.ffMeasure;
         propList.measure = obj.measure;
         propList.absTol = obj.absTol;
         propList.relTol = obj.relTol;
         
         % Willy added (display following)
         %          propList.mmin=obj.mmin;
         %          propList.mmax=obj.mmax;
         
         if obj.nInit ~= obj.def_nInit
            propList.nInit = obj.nInit;
         end
         if obj.nMax ~= obj.def_nMax
            propList.nMax = obj.nMax;
         end
         if obj.nMu ~= obj.def_nMu
            propList.nMu = obj.nMu;
         end
         if obj.nf ~= obj.def_nf
            propList.nf= obj.nf;
         end
         if numel(obj.trueMuCV)
            propList.trueMuCV = obj.trueMuCV;
         end
      end
      
   end
   
end

