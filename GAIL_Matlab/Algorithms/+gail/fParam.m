classdef fParam < gail.errorParam
   %GAIL.FPARAM is a class containing the parameters related to
   %algorithms that act on functions of x
   %   This class contains the function, its domain, etc.
   
   properties
      f %function defined on some domain, often an interval or box
      domain %domain one or several dimensions
      nInit %initial sample size
      nMax %maximum sample size
      
   end
   
    properties (Dependent = true)
       d %number of variables
       nfOut %number of outputs of f
    end
   
   properties (Hidden, SetAccess = private)
      def_f = @(x) sum(x.^2,2) %default function
      def_domain = [0; 1]; %default domain
      def_nInit = 100 %default initial number of samples
      def_nMax = 1e7 %default maximum sample size
   end
   
   
   methods
      
      % Creating a meanYParam process
      function obj = fParam(varargin)
         %this constructor essentially parses inputs
         %the parser will look for the following in order
         %  # a copy of a fParam object
         %  # a function
         %  # a domain
         %  # a structure
         %  # numbers: absTol, relTol,
         %  # name-value pairs
         
         start = 1;
         objInp = 0;
         fInp = 0;
         domainInp = 0;
         structInp = 0;
         if nargin %there are inputs to parse and assign
            if isa(varargin{start},'gail.fParam') 
               %the first input is a meanYParam object so copy it
               objInp = start;
               start = start + 1;
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
                     end
                  end
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
            obj.f = val.f; %copy random number generator
            obj.domain = val.domain; %copy uncertainty
            obj.nInit = val.nInit; %copy sample size for sigma
            obj.nMax = val.nMax; %copy maximum sample size
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
         f_addParamVal(p,'f',obj.def_f);
         f_addParamVal(p,'domain',obj.def_domain);
         f_addParamVal(p,'nInit',obj.def_nInit);
         f_addParamVal(p,'nMax',obj.def_nMax);
         
         if structInp
            parse(p,varargin{start:end},varargin{structInp}) 
            %parse inputs with a structure
         else
            parse(p,varargin{start:end}) %parse inputs w/o structure
         end
         struct_val = p.Results; %store parse inputs as a structure
         
         %Assign values of structure to corresponding class properties
         if fInp
            obj.f = varargin{fInp}; %assign function
         else
            obj.f = struct_val.f;
         end
         if domainInp
            obj.domain = varargin{domainInp}; %assign function
         else
            obj.domain = struct_val.domain;
         end
         obj.nInit = struct_val.nInit;
         obj.nMax = struct_val.nMax;

         
      end %of constructor
     
      function set.f(obj,val)
         validateattributes(val, {'function_handle'}, {})
         obj.f = val;
      end
      
      function set.domain(obj,val)
         validateattributes(val, {'numeric'}, {'2d'})
         assert(numel(val) > 1)
         if any(size(val) == 1) %domain is an interval
            val = val(:)';
         end
         obj.domain = val;
      end
                       
      function set.nInit(obj,val)
         validateattributes(val, {'numeric'}, {'scalar','positive','integer'})
         obj.nInit = val;
      end
                             
      function set.nMax(obj,val)
         validateattributes(val, {'numeric'}, {'scalar','positive','integer'})
         obj.nMax = val;
      end
                                                                                
      function val = get.d(obj)
         val = size(obj.domain,2); 
      end         
                       
      function val = get.nfOut(obj)
         val = numel(obj.f(obj.domain(1,:))); 
      end      
       
                       
   end
   

   
   
end

