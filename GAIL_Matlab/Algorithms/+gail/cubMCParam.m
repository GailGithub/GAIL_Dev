classdef cubMCParam<gail.cubParam
   properties
   Y %function to integrate
   alpha %uncertainty
   nSig %sample size to estimate variance
   nYOut %number of Y for each mean
   nY
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
   whichParse = [objInp structInp start:nargin];
   whichParse = whichParse(whichParse > 0);
   obj@gail.cubParam(varargin{whichParse});
   if objInp
            val = varargin{objInp}; %first input
            obj.Y = val.Y; %copy random number generator
            obj.alpha = val.alpha; %copy uncertainty
            obj.nSig = val.nSig; %copy sample size for sigma
            obj.nY =obj.nf; %copy number of Y values per mu
            obj.nYOut=obj.fun.nfOut;
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
           start = start + 2; % ???  to account for the two tolerances already parsed
         end
         f_addParamVal(p,'Y',obj.def_Y);
         f_addParamVal(p,'alpha',obj.def_alpha);
         f_addParamVal(p,'nSig',obj.def_nSig);
         f_addParamVal(p,'nY',obj.def_nY);
         
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
         if isfield(struct_val,'nY')
            obj.nY = struct_val.nY;
         end
         
   end %end of constructor
   
    % f = @(x) prod(x,2);
    % q = meanMC_CLT(@(n)f(rand(n,2)),1e-4,0); exactsol = 1/4; 
   function val=get.Y(obj)
      if strcmp(obj.measureType,'uniform')
         val = @(n)obj.ff(rand(n,obj.fun.d));

      elseif strcmp(obj.measureType,'normal') || strcmp(obj.measureType,'Gaussian')
         val = @(n)obj.ff(rand(n,obj.fun.d));
      end
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
      
       function set.nY(obj,val)
         validateattributes(val, {'numeric'}, {'positive','integer'})
         obj.nY = val;
       end
       
   end
end