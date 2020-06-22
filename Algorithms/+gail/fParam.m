classdef fParam < handle & matlab.mixin.CustomDisplay
   %GAIL.FPARAM is a class containing the parameters related to
   %algorithms that act on functions of x
   %   This class contains the function, its domain, etc.
   %
   % Example 1.  Default values
   % >> fParamObj = gail.fParam
   % fParamObj = ***
   %
   %              f: @(x)sum(x.^2,2)
   %         domain: [2***1 double]
   %
   %
   % Example 2. Function on 2-D unit cube
   % >> fParamObj = gail.fParam(@(x) sum(x.^3.2),[0 0; 1 1])
   % fParamObj = ***
   %
   %              f: @(x)sum(x.^3.2)
   %         domain: [2***2 double]
   %
   %
   % Example 3. Ball domain
   % >> fParamObj = gail.fParam(@(x) sum(x.^3.2),[0 0; 1 1],'ball')
   % fParamObj = ***
   %
   %              f: @(x)sum(x.^3.2)
   %         domain: [2***2 double]
   %     domainType: 'ball'
   %
   %
   % Example 4. Using name/value pairs
   % fParamObj = gail.fParam('domain', [0 0; 1 1], 'f', @(x) sum(x.^3.2), 'relTol', 0.1)
   % fParamObj = ***
   %
   %              f: @(x)sum(x.^3.2)
   %         domain: [2***2 double]
   %
   %
   % Example 5. Using structure for input
   % >> inpStruct.f = @(x) sin(sum(x,2));
   % >> inpStruct.domain = [zeros(1,4); ones(1,4)];
   % >> inpStruct.relTol = 0.1;
   % >> fParamObj = gail.fParam(inpStruct)
   % fParamObj = ***
   %
   %              f: @(x)sin(sum(x,2))
   %         domain: [2***4 double]
   %
   %
   % Example 6. Using structure for input and numbers, structure takes
   % precedence
   % >> inpStruct.f = @(x) sin(sum(x,2));
   % >> inpStruct.domain = [zeros(1,4); ones(1,4)];
   % >> inpStruct.relTol = 0.1;
   % >> fParamObj = gail.fParam(inpStruct,0.0001,0.01)
   % fParamObj = ***
   %
   %              f: @(x)sin(sum(x,2))
   %         domain: [2***4 double]
   %
   %
   % Example 7. Coping an fParam object and changing properties
   % >> newfParamObj = gail.fParam(fParamObj,'domainType','sphere')
   % newfParamObj = ***
   %
   %              f: @(x)sin(sum(x,2))
   %         domain: [2***4 double]
   %     domainType: 'sphere'
   %
   %
   % Author:  Fred J. Hickernell

   properties
      f %function defined on some domain, often an interval or box
      domain %domain one or several dimensions
      domainType %domain type
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
      def_nInit = 1e3 %default initial number of samples
      def_nMax = 1.5e6 %default maximum sample size
      allowedDomains = {'box', ... %a hyperbox
         'ball', ... %solid ball
         'sphere',...
         'ball-from-normal', ...
         'ball-from-cube', ...
         'sphere-from-normal',...
         'sphere-from-cube'
         } %hollow sphere
   end


   methods

      % Creating a fParam process
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

         start = 1; %index to begin to parse
         useDefaults = true; %true unless copying an fParam object, then false
         objInp = 0; %where is the an object in the class
         fInp = 0; %where is the input function
         domainInp = 0; %where is the domain
         domainTypeInp = 0; %where is the domain type
         structInp = 0; %where is the structure
         if nargin %there are inputs to parse and assign
            if isa(varargin{start},'gail.fParam')
               %the first input is a fParam object so copy it
               objInp = start;
               start = start + 1;
            end
            if nargin >= start
               if isstruct(varargin{start}) %next input is a structure containing f
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

         done = false; %not finished parsing
         if objInp
            val = varargin{objInp}; %first input
            obj.f = val.f; %copy function
            obj.domain = val.domain; %copy domain
            obj.domainType = val.domainType; %copy domain type
            obj.nInit = val.nInit; %copy initial sample size
            obj.nMax = val.nMax; %copy maximum sample size
            useDefaults = false;
         end

         %Now begin to parse inputs
         p = inputParser; %construct an inputParser object
         p.KeepUnmatched = true; %ignore those that do not match
         p.PartialMatching = false; %don'try a partial match
         p.StructExpand = true; %expand structures
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
         if ~useDefaults %remove defaults if copying fParam object
            struct_val = rmfield(struct_val,p.UsingDefaults);
         end

         %Assign values of structure to corresponding class properties
         if fInp
            obj.f = varargin{fInp}; %assign function
         elseif isfield(struct_val,'f')
            if any(strcmp(p.UsingDefaults,'f'))
               warning('GAIL:fParam:noFunctionInput','No function input, default used.')
            end
            obj.f = struct_val.f;
         end
         if domainInp
            obj.domain = varargin{domainInp}; %assign function
         elseif isfield(struct_val,'domain')
            obj.domain = struct_val.domain;
         end
         if domainTypeInp
            obj.domainType = varargin{domainTypeInp}; %assign function
         elseif isfield(struct_val,'domainType')
            obj.domainType = struct_val.domainType;
         end
         if isfield(struct_val,'nInit')
            obj.nInit = struct_val.nInit;
         end
         if isfield(struct_val,'nMax')
            obj.nMax = struct_val.nMax;
         end

      end %of constructor

      function set.f(obj,val)
         validateattributes(val, {'function_handle'}, {'nonempty'})
         obj.f = val;
      end

      function set.domain(obj,val)
         validateattributes(val, {'numeric'}, {'2d'})
         assert(numel(val) > 1)
         if size(val, 1) == 1 %domain is an interval
            val = val(:);
         end
         obj.domain = val(1:2,:);
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
         % val = numel(obj.f(obj.domain(1,:)));
      end

   end

  methods (Access = protected)
      function checkDomainType(obj,inval)
         assert(any(strcmp(inval,obj.allowedDomains)))
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
         propList = struct('f', obj.f, ...
            'domain', obj.domain);
         if ~strcmp(obj.domainType,obj.def_domainType)
            propList.domainType = obj.domainType;
         end
         if obj.nInit ~= obj.def_nInit
            propList.nInit = obj.nInit;
         end
         if obj.nMax ~= obj.def_nMax
            propList.nMax = obj.nMax;
         end
      end
  end

end
