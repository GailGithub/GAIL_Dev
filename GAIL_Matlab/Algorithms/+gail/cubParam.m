classdef cubParam < gail.fParam
   %GAIL.CUBPARAM is a class containing the parameters related to
   %algorithms that find the mean of a random variable
   %   This class contains the number of integrands with the same integral,
   %   etc.
   
   properties
      measure %measure against which to integrate
      inflate %inflation factor for bounding the error
      nMu %number of integrals for solution function
      nf %number of f for each integral
      trueMuCV %true integral for control variates
   end
   
    properties (Dependent = true)
       ff %function after variable transformation
       nCV %number of control variates
       volume %volume of the domain with respect to the measure
    end
   
   properties (Hidden, SetAccess = private)
      def_measure = 'uniform'
      def_inflate = @(m) 5 * 2^-m %default inflation factor
      def_nMu = 1 %default number of integrals
      def_nf = 1 %default number of Y per integral
      def_trueMuCV = [] %default true integrals for control variates
      allowedMeasures = {'uniform', ... %over a hyperbox volume of domain is one
         'Lebesgue', ... %like uniform, but integral over domain is the volume of the domain
         'Gaussian', 'normal', ... %these are the same
         'Uniform ball', ... %analogous 
         'Lebesgue ball', ... %to hyperbox
         'Uniform sphere', ... %analogous 
         'Lebesgue sphere', ... %to hyperbox
         }
   end
   
   
   methods
      
      % Creating a cubParam process
      function obj = cubParam(varargin)
         %this constructor essentially parses inputs
         %the parser will look for the following in order
         %  # a copy of a cubParam object
         %  # a function
         %  # a domain
         %  # a measure
         %  # a structure
         %  # numbers: absTol, relTol,
         %  # name-value pairs
         
         start = 1;
         objInp = 0;
         fInp = 0;
         domainInp = 0;
         measureInp = 0;
         structInp = 0;
         if nargin %there are inputs to parse and assign
            if isa(varargin{start},'gail.cubParam') 
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
                        if nargin >= start
                           if ischar(varargin{start})
                              measureInp = start;
                              start = start + 1;
                           end
                        end
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
         whichParse = [objInp fInp domainInp structInp start:nargin];
         whichParse = whichParse(whichParse > 0);
         obj@gail.fParam(varargin{whichParse});

         if objInp
            val = varargin{objInp}; %first input
            obj.measure = val.measure; %copy integration measure
            obj.inflate = val.inflate; %copy inflation factor
            obj.nMu = val.nMu; %copy the number of integrals
            obj.nf = val.nf; %copy number of functions for each integral
            obj.trueMuCV = val.trueMuCV; %copy true means of control variates
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
         f_addParamVal(p,'measure',obj.def_measure);
         f_addParamVal(p,'inflate',obj.def_inflate);
         f_addParamVal(p,'nMu',obj.def_nMu);
         f_addParamVal(p,'nf',obj.def_nf);
         f_addParamVal(p,'trueMuCV',obj.def_trueMuCV);
         
         if structInp
            parse(p,varargin{start:end},varargin{structInp}) 
            %parse inputs with a structure
         else
            parse(p,varargin{start:end}) %parse inputs w/o structure
         end
         struct_val = p.Results; %store parse inputs as a structure
         
         %Assign values of structure to corresponding class properties
         if measureInp
            obj.measure = varargin{measureInp}; %assign measure
         else
            obj.measure = struct_val.measure;
         end
         obj.inflate = struct_val.inflate;
         obj.nMu = struct_val.nMu;
         obj.nf = struct_val.nf;
         obj.trueMuCV = struct_val.trueMuCV;
       
      end %of constructor
     
      function set.measure(obj,val)
         validateattributes(val, {'char'}, {})
         obj.measure = checkMeasure(obj,val);
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
      
      function val = get.nCV(obj)
         val = obj.nfOut - sum(obj.nf); 
      end        
      
      function val = get.volume(obj) %volume of the domain
         if any(strcmp(obj.measure,{'uniform', 'normal'}))
            val = 1;
         elseif strcmp(obj.measure, {'Lebesgue'})
            val = prod(diff(obj.domain,1),2);
         end
      end
            
      function val = get.ff(obj)
         if strcmp(obj.measure,'uniform')
            val = obj.f;
         elseif strcmp(obj.measure,'Lebesgue')
            val = @(x) obj.f(bsxfun(@times, diff(obj.domain,1), ...
               bsxfun(@minus, x, obj.domain(1,:))));
         elseif strcmp(obj.measure, 'normal')
            val = @(x) obj.f(gail.stdnorminv(x));
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
            assert(all(obj.domain(1,:) == '-Inf') && ...
               all(obj.domain(2,:) == 'Inf'))
         end
         if any(strcmp(outval,{'uniform','Lebesgue'})) %domain must be finite
            assert(all(all(isfinite(obj.domain))))
         end
      end
   
      function outval = setTrueMuCVDim(obj,inval)
         assert(numel(inval) == obj.nCV)
         outval = inval(:)';
      end
      
   end

   
   
end

