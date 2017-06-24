classdef fParam < gail.errorParam
   %GAIL.FPARAM is a class containing the parameters related to
   %algorithms that act on functions of x
   %   This class contains the function, its domain, etc.
   %
   % Example 1.  Default values
   % >> fParamObj = gail.fParam 
   % fParamObj = 
   %   fParam with properties:
   % 
   %              f: @(x)sum(x.^2,2)
   %         domain: [2×1 double]
   %     domainType: 'box'
   %          nInit: 100
   %           nMax: 10000000
   %              d: 1
   %          nfOut: 1
   %         absTol: 0.010000000000000
   %         relTol: 0
   %         solFun: @(mu)mu
   %       solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
   %
   %
   % Example 2. Function on 2-D unit cube
   % >> fParamObj = gail.fParam(@(x) sum(x.^3.2),[0 0; 1 1])
   % fParamObj = 
   %   fParam with properties:
   % 
   %              f: @(x)sum(x.^3.2)
   %         domain: [2×2 double]
   %     domainType: 'box'
   %          nInit: 100
   %           nMax: 10000000
   %              d: 2
   %          nfOut: 1
   %         absTol: 0.010000000000000
   %         relTol: 0
   %         solFun: @(mu)mu
   %       solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
   % 
   %
   % Example 3. Ball domain
   % >> fParamObj = gail.fParam(@(x) sum(x.^3.2),[0 0; 1 1],'ball')
   % fParamObj = 
   %   fParam with properties:
   % 
   %              f: @(x)sum(x.^3.2)
   %         domain: [2×2 double]
   %     domainType: 'ball'
   %          nInit: 100
   %           nMax: 10000000
   %              d: 2
   %          nfOut: 1
   %         absTol: 0.010000000000000
   %         relTol: 0
   %         solFun: @(mu)mu
   %       solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
   %
   %
   % Example 4. Using name/value pairs
   % fParamObj = gail.fParam('domain', [0 0; 1 1], 'f', @(x) sum(x.^3.2), 'relTol', 0.1)
   % fParamObj = 
   %   fParam with properties:
   % 
   %              f: @(x)sum(x.^3.2)
   %         domain: [2×2 double]
   %     domainType: 'box'
   %          nInit: 100
   %           nMax: 10000000
   %              d: 2
   %          nfOut: 1
   %         absTol: 0.010000000000000
   %         relTol: 0.100000000000000
   %         solFun: @(mu)mu
   %       solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
   % 
   %
   % Example 5. Using structure for input
   % >> inpStruct.f = @(x) sin(sum(x,2));
   % >> inpStruct.domain = [zeros(1,4); ones(1,4)];
   % >> inpStruct.relTol = 0.1;
   % >> fParamObj = gail.fParam(inpStruct)
   % fParamObj = 
   %   fParam with properties:
   % 
   %              f: @(x)sin(sum(x,2))
   %         domain: [2×4 double]
   %     domainType: 'box'
   %          nInit: 100
   %           nMax: 10000000
   %              d: 4
   %          nfOut: 1
   %         absTol: 0.010000000000000
   %         relTol: 0.100000000000000
   %         solFun: @(mu)mu
   %       solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
   %
   %   
   % Example 6. Using structure for input and numbers, structure takes
   % precedence
   % >> inpStruct.f = @(x) sin(sum(x,2));
   % >> inpStruct.domain = [zeros(1,4); ones(1,4)];
   % >> inpStruct.relTol = 0.1;
   % >> fParamObj = gail.fParam(inpStruct,0.0001,0.01)
   % fParamObj = 
   %   fParam with properties:
   % 
   %              f: @(x)sin(sum(x,2))
   %         domain: [2×4 double]
   %     domainType: 'box'
   %          nInit: 100
   %           nMax: 10000000
   %              d: 4
   %          nfOut: 1
   %         absTol: 1.000000000000000e-04
   %         relTol: 0.100000000000000
   %         solFun: @(mu)mu
   %       solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
   %
   %   
   % Author:  Fred J. Hickernell
   
   
   
   properties
      f %function defined on some domain, often an interval or box
      domain %domain one or several dimensions
      domainType
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
      def_domainType = 'box'; %default domain type
      def_nInit = 100 %default initial number of samples
      def_nMax = 1e7 %default maximum sample size
      allowedDomains = {'box', ... %a hyperbox
         'ball', ... %solid ball
         'sphere'} %hollow sphere
   end
   
   
   methods
      
      % Creating a meanYParam process
      function obj = fParam(varargin)
         %this constructor essentially parses inputs
         %the parser will look for the following in order
         %  # a copy of a fParam object
         %  # a structure
         %  # a function
         %  # a domain
         %  # a domain Type
         %  # numbers: absTol, relTol
         %  # name-value pairs
         
         start = 1;
         objInp = 0;
         fInp = 0;
         domainInp = 0;
         domainTypeInp = 0;
         structInp = 0;
         if nargin %there are inputs to parse and assign
            if isa(varargin{start},'gail.fParam') 
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
                        start = start + 1;
                        if nargin >= start 
                           if ischar(varargin{start}) %next input is domain type
                              domainTypeInp = start;
                              start = start+ 1;
                           end
                        end
                     end
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
            obj.f = val.f; %copy function
            obj.domain = val.domain; %copy domain
            obj.domainType = val.domainType; %copy domain type
            obj.nInit = val.nInit; %copy initial sample size
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
               parseRange = start:nargin;
               done = true;
            end
         end
         if ~done %then nothingleft or just numbers
           f_addParamVal = @addOptional;
           parseRange = []; %nothing to parse here if just numbers
         end
         if ~fInp
            f_addParamVal(p,'f',obj.def_f);
         end
         if ~domainInp
            f_addParamVal(p,'domain',obj.def_domain);
         end
         if ~domainTypeInp
            f_addParamVal(p,'domainType',obj.def_domainType);
         end
         f_addParamVal(p,'nInit',obj.def_nInit);
         f_addParamVal(p,'nMax',obj.def_nMax);
         
         if structInp
            parse(p,varargin{parseRange},varargin{structInp}) 
            %parse inputs with a structure
         else
            parse(p,varargin{parseRange}) %parse inputs w/o structure
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
         if domainTypeInp
            obj.domainType = varargin{domainTypeInp}; %assign function
         else
            obj.domainType = struct_val.domainType;
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
            val = val(:);
         end
         obj.domain = val;
      end
                       
      function set.domainType(obj,val)
         validateattributes(val, {'char'}, {})
         checkDomainType(obj,val)
         obj.domainType = val;
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
   
  methods (Access = protected)
   
      function checkDomainType(obj,inval)
         assert(any(strcmp(inval,obj.allowedDomains)))
      end
         
   end

   
   
end

