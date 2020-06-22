classdef cubMCParam < gail.cubParam
   properties
      alpha %uncertainty
      nSig %sample size to estimate variance
   end

   properties (Dependent)
      Y %random variable to integrate
      nY %number of Y output
      nYOut %number of Y for each mean
   end

   properties (Hidden, SetAccess = private)
      def_Y = @(n) rand(n,1) %default random number generator
      def_alpha = 0.01 %default uncertainty
      def_nSig = 1000 %default uncertainty
      def_nYOut = 1 %default number of Y per mean
      def_nY=1 %default number of Y
   end
   
   methods
   function obj = cubMCParam(varargin)
      %this constructor essentially parses inputs
      %the parser will look for the following in order
      %  # a copy of a meanYParam object
      %  # a function
      %  # a structure
      %  # numbers: absTol, relTol, alpha,
      %  # name-value pairs
      start = 1;
      objInp = 0;
      if nargin %there are inputs to parse and assign
         if isa(varargin{start},'gail.cubMCParam') 
            %the first input is a meanYParam object so copy it
            objInp = start;
            start = start + 1;
         end
      end
      whichParse = [objInp start:nargin];
      whichParse = whichParse(whichParse > 0);
      obj@gail.cubParam(varargin{whichParse});
      if objInp
         val = varargin{objInp}; %first input
%         obj.Y = val.Y; %copy random number generator
         obj.alpha = val.alpha; %copy uncertainty
         obj.nSig = val.nSig; %copy sample size for sigma
%          obj.nY = val.nY; %copy number of Y values per mu
%          obj.nYOut = val.nYOut;
         return
      end

      %Now begin to parse inputs
      p = inputParser; %construct an inputParser object
      p.KeepUnmatched = true; %ignore those that do not match
      p.PartialMatching = false; %don'try a partial match
      p.StructExpand = true; %expand structures
      f_addParamVal = @addOptional;
      parseRange = []; %nothing to parse here if just numbers
      if nargin >= start
         if ischar(varargin{start})
            %there may be input string/value pairs or a structure
            MATLABVERSION = gail.matlab_version;
            if MATLABVERSION >= 8.3
               f_addParamVal = @addParameter;
            else
               f_addParamVal = @addParamValue;
            end
            parseRange = 1:nargin;
         end
      end
      f_addParamVal(p,'alpha',obj.def_alpha);
      f_addParamVal(p,'nSig',obj.def_nSig);
      parse(p,varargin{parseRange}) %parse inputs w/o structure
      struct_val = p.Results; %store parse inputs as a structure

      %Assign values of structure to corresponding class properties
      if isfield(struct_val,'alpha')
         obj.alpha = struct_val.alpha;
      end
      if isfield(struct_val,'nSig')
         obj.nSig = struct_val.nSig;
      end

   end %end of constructor
   
    % f = @(x) prod(x,2);
    % q = meanMC_CLT(@(n)f(rand(n,2)),1e-4,0); exactsol = 1/4; 
   function set.alpha(obj,val)
      validateattributes(val, {'numeric'}, {'scalar','nonnegative', ...
         '<', 1})
      obj.alpha = val;
   end

   function set.nSig(obj,val)
      validateattributes(val, {'numeric'}, {'scalar','positive','integer'})
      obj.nSig = val;
   end

   function val = get.Y(obj)
      if strcmp(obj.measureType,'uniform')
         val = @(n) obj.ff(rand(n,obj.fun.d));
      elseif strcmp(obj.measureType,'normal') || strcmp(obj.measureType,'Gaussian')
         val = @(n) obj.ff(randn(n,obj.fun.d));
      end
   end
   
   function val = get.nY(obj)
      val = obj.nf;
   end
       
   function val = get.nYOut(obj)
      val = obj.fun.nfOut;
   end
       
   end
end