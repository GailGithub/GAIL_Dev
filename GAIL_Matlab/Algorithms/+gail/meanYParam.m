classdef meanYParam < gail.errorParam
   %GAIL.MEANYPARAM is a class containing the parameters related to
   %algorithms that find the mean of a random variable
   %   This class contains the random number generator, the uncertainty
   %
   % Example 1. Construct a meanYParam object with default parameters
   % >> meanYParamObj = gail.meanYParam
   % meanYParamObj = 
   %   meanYParam with properties:
   % 
   %            Y: @(n)rand(n,1)
   %        alpha: 0.010000000000000
   %         nSig: 1000
   %      inflate: 1.200000000000000
   %         nMax: 100000000
   %          nMu: 1
   %           nY: 1
   %     trueMuCV: [1×0 double]
   %          nCV: 0
   %        nYOut: 1
   %       absTol: 0.010000000000000
   %       relTol: 0
   %       solFun: @(mu)mu
   %     solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
   %
   %
   % Example 2. Construct a cubParam object with properly ordered inputs
   % >> meanYParamObj = gail.meanYParam(@(n) sum(rand(n,4).^3,2),0.001)
   % meanYParamObj = 
   %   meanYParam with properties:
   % 
   %            Y: @(n)sum(rand(n,4).^3,2)
   %        alpha: 0.010000000000000
   %         nSig: 1000
   %      inflate: 1.200000000000000
   %         nMax: 100000000
   %          nMu: 1
   %           nY: 1
   %     trueMuCV: [1×0 double]
   %          nCV: 0
   %        nYOut: 1
   %       absTol: 1.000000000000000e-03
   %       relTol: 0
   %       solFun: @(mu)mu
   %     solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
   %
   %
   % Example 3. Using name/value pairs
   % >> meanYParamObj = gail.meanYParam('nSig', 1e4, 'Y', @(x) sin(sum(x.^3,2)), 'relTol', 0.1)
   % meanYParamObj = 
   %   meanYParam with properties:
   % 
   %            Y: @(x)sin(sum(x.^3,2))
   %        alpha: 0.010000000000000
   %         nSig: 10000
   %      inflate: 1.200000000000000
   %         nMax: 100000000
   %          nMu: 1
   %           nY: 1
   %     trueMuCV: [1×0 double]
   %          nCV: 0
   %        nYOut: 1
   %       absTol: 0.010000000000000
   %       relTol: 0.100000000000000
   %       solFun: @(mu)mu
   %     solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
   %
   %
   % Example 4. Using a structure for input
   % >> inpStruct.Y = @(x) sin(sum(x,2));
   % >> inpStruct.nSig = 1e4;
   % >> inpStruct.relTol = 0.1;
   % >> meanYParamObj = gail.meanYParam(inpStruct)
   % meanYParamObj = 
   %   meanYParam with properties:
   % 
   %            Y: @(x)sin(sum(x,2))
   %        alpha: 0.010000000000000
   %         nSig: 10000
   %      inflate: 1.200000000000000
   %         nMax: 100000000
   %          nMu: 1
   %           nY: 1
   %     trueMuCV: [1×0 double]
   %          nCV: 0
   %        nYOut: 1
   %       absTol: 0.010000000000000
   %       relTol: 0.100000000000000
   %       solFun: @(mu)mu
   %     solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
   %
   %
   % Example 5. Copying a meanYParam object and changing some properties
   % >> NewMeanYParamObj = gail.meanYParam(meanYParamObj,'Y',@(n) rand(n,3))
   % NewMeanYParamObj = 
   %   meanYParam with properties:
   % 
   %            Y: @(n)rand(n,3)
   %        alpha: 0.010000000000000
   %         nSig: 10000
   %      inflate: 1.200000000000000
   %         nMax: 100000000
   %          nMu: 1
   %           nY: 1
   %     trueMuCV: [1×0 double]
   %          nCV: 2
   %        nYOut: 3
   %       absTol: 0.010000000000000
   %       relTol: 0.100000000000000
   %       solFun: @(mu)mu
   %     solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
   %
   %
   % Author: Fred J. Hickernell

   properties
      Y %random number generator
      alpha %uncertainty
      nSig %sample size to estimate variance
      inflate %inflation factor for bounding standard deviation
      nMax %maximum sample size allowed
      nMu %number of means for solution function
      nY %number of Y for each mean
      trueMuCV %true mean of control variates
   end
   
    properties (Dependent = true)
       nCV %number of control variates
       nYOut %number of Y outputs
    end
   
   properties (Hidden, SetAccess = private)
      def_Y = @(n) rand(n,1) %default random number generator
      def_alpha = 0.01 %default uncertainty
      def_nSig = 1000 %default uncertainty
      def_inflate = 1.2 %default inflation factor
      def_nMax = 1e8 %default maximum sample size
      def_nMu = 1 %default number of means
      def_nY = 1 %default number of Y per mean
      def_trueMuCV = [] %default true means of control variates
   end
   
   
   methods
      
      % Creating a meanYParam process
      function obj = meanYParam(varargin)
         %this constructor essentially parses inputs
         %the parser will look for the following in order
         %  # a copy of a meanYParam object
         %  # a function
         %  # a structure
         %  # numbers: absTol, relTol, alpha,
         %  # name-value pairs
         
         start = 1;
         useDefaults = true;
         objInp = 0;
         YInp = 0;
         structInp = 0;
         if nargin %there are inputs to parse and assign
            if isa(varargin{start},'gail.meanYParam') 
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
               if gail.isfcn(varargin{start}) %next input is the function Y
                  YInp = start;
                  start = start + 1;
               end
            end
         end
         
         %Parse errorParam properties
         whichParse = [objInp structInp start:nargin];
         whichParse = whichParse(whichParse > 0);
         obj@gail.errorParam(varargin{whichParse});

         if objInp
            val = varargin{objInp}; %first input
            obj.Y = val.Y; %copy random number generator
            obj.alpha = val.alpha; %copy uncertainty
            obj.nSig = val.nSig; %copy sample size for sigma
            obj.inflate = val.inflate; %copy inflation factor
            obj.nMax = val.nMax; %copy maximum sample size
            obj.nMu = val.nMu; %copy number of mu values
            obj.nY = val.nY; %copy number of Y values per mu
            obj.trueMuCV = val.trueMuCV; %copy true means of control variates
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
               done = true;
            end
         end
         if ~done %then nothingleft or just numbers
           f_addParamVal = @addOptional;
           start = start + 2; %to account for the two tolerances already parsed
         end
         f_addParamVal(p,'Y',obj.def_Y);
         f_addParamVal(p,'alpha',obj.def_alpha);
         f_addParamVal(p,'nSig',obj.def_nSig);
         f_addParamVal(p,'inflate',obj.def_inflate);
         f_addParamVal(p,'nMax',obj.def_nMax);
         f_addParamVal(p,'nMu',obj.def_nMu);
         f_addParamVal(p,'nY',obj.def_nY);
         f_addParamVal(p,'trueMuCV',obj.def_trueMuCV);
         
         if structInp
            parse(p,varargin{start:end},varargin{structInp}) 
            %parse inputs with a structure
         else
            parse(p,varargin{start:end}) %parse inputs w/o structure
         end
         struct_val = p.Results; %store parse inputs as a structure
         if ~useDefaults %remove defaults if copying cubParam object
            struct_val = rmfield(struct_val,p.UsingDefaults);
         end
         
         %Assign values of structure to corresponding class properties
         if YInp
            obj.Y = varargin{YInp}; %assign function
         elseif isfield(struct_val,'Y')
            obj.Y = struct_val.Y;
         end
         if isfield(struct_val,'alpha')
            obj.alpha = struct_val.alpha;
         end
         if isfield(struct_val,'nSig')
            obj.nSig = struct_val.nSig;
         end
         if isfield(struct_val,'inflate')
            obj.inflate = struct_val.inflate;
         end
         if isfield(struct_val,'nMax')
            obj.nMax = struct_val.nMax;
         end
         if isfield(struct_val,'nMu')
            obj.nMu = struct_val.nMu;
         end
         if isfield(struct_val,'nY')
            obj.nY = struct_val.nY;
         end
         if isfield(struct_val,'trueMuCV')
            obj.trueMuCV = struct_val.trueMuCV;
         end

         
      end %of constructor
     
      function set.Y(obj,val)
         validateattributes(val, {'function_handle'}, {})
         obj.Y = val;
      end
      
      function set.alpha(obj,val)
         validateattributes(val, {'numeric'}, {'scalar','nonnegative', ...
            '<', 1})
         obj.alpha = val;
      end
                       
      function set.nSig(obj,val)
         validateattributes(val, {'numeric'}, {'scalar','positive','integer'})
         obj.nSig = val;
      end
                       
      function set.inflate(obj,val)
         validateattributes(val, {'numeric'}, {'scalar','>',1})
         obj.inflate = val;
      end
      
      function set.nMax(obj,val)
         validateattributes(val, {'numeric'}, {'scalar','positive','integer'})
         obj.nMax = val;
      end
      
      function set.nMu(obj,val)
         validateattributes(val, {'numeric'}, {})
         obj.nMu = val;
      end
                       
      function set.nY(obj,val)
         validateattributes(val, {'numeric'}, {'positive','integer'})
         obj.nY = val;
      end
                       
      function set.trueMuCV(obj,val)
         validateattributes(val, {'numeric'}, {})
         obj.trueMuCV = setTrueMuCVDim(obj,val);
      end
      
      function val = get.nCV(obj)
         val = obj.nYOut - sum(obj.nY); 
      end         
                       
      function val = get.nYOut(obj)
         val = numel(obj.Y(1)); 
      end         
                       
     
      
   end
   
   methods (Access = protected)
   
      function outval = setTrueMuCVDim(obj,inval)
         assert(numel(inval) == obj.nCV)
         outval = inval(:)';
      end
   
   end

   
   
end

