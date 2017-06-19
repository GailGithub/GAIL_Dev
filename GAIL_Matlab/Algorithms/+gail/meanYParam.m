classdef meanYParam < gail.errorParam
   %GAIL.MEANYPARAM is a class containing the parameters related to
   %algorithms that find the mean of a random variable
   %   This class contains the random number generator, the uncertainty
   
   properties
      Y %random number generator
      alpha %uncertainty
      nSig %sample size to estimate variance
      inflate %inflation factor for bounding standard deviation
      nMax %maximum sample size allowed
   end
   
   properties (Hidden, SetAccess = private)
      def_Y = @(n) rand(n,1) %default random number generator
      def_alpha = 0.01 %default uncertainty
      def_nSig = 1000 %default uncertainty
      def_inflate = 1.2 %default inflation factor
      def_nMax = 1e8 %default maximum sample size
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
               if gail.isfcn(varargin{start}) %next input is the function Y
                  YInp = start;
                  start = start + 1;
               end
               if nargin >= start
                  if isstruct(varargin{start}) %next input is a structure containing Y
                     structInp = start;
                     start = start + 1;
                  end
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
            obj.nMax = val.nMax; %copy maximum ample size
            return
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
         if structInp
            parse(p,varargin{start:end},varargin{structInp}) 
            %parse inputs with a structure
         else
            parse(p,varargin{start:end}) %parse inputs w/o structure
         end
         struct_val = p.Results; %store parse inputs as a structure
         
         %Assign values of structure to corresponding class properties
         obj.alpha = struct_val.alpha;
         obj.nSig = struct_val.nSig;
         obj.inflate = struct_val.inflate;
         obj.nMax = struct_val.nMax;
         if YInp
            obj.Y = varargin{YInp}; %assign function
         else
            obj.Y = struct_val.Y;
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
                       
     
      
   end
   
end

