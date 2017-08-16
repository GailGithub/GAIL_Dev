classdef cubMeanParam < handle & matlab.mixin.CustomDisplay
   %GAIL.CUBMEANPARAM is a class containing the parameters related to
   %algorithms that act on functions of x
   %   This class contains the function, its domain, etc.
   %
   % Example 1. Construct a cubMeanParam object with default parameters
   % >> cubMeanParamObj = gail.cubMeanParam
   % cubMeanParamObj = 
   %   cubMeanParam with properties:
   % 
   %     nInit: 1024
   %      nMax: 16777216
   %
   %
   % Example 2. Using name/value pairs
   % cubMeanParamObj = gail.cubMeanParam('nInit', 100)
   % cubMeanParamObj = 
   %   cubMeanParam with properties:
   % 
   %     nInit: 100
   %      nMax: 16777216
   % 
   %
   % Example 2. Using structure for input
   % >> inpStruct.nInit = 200;
   % >> cubMeanParamObj = gail.cubMeanParam(inpStruct)
   % cubMeanParamObj = 
   %   cubMeanParam with properties:
   % 
   %     nInit: 200
   %      nMax: 16777216
   %
   %   
   % Example 3. Copying an cubMeanParamObj object and changing properties
   % >> newcubMeanParamObj = gail.cubMeanParam(cubMeanParamObj,'nInit',300)
   % newcubMeanParamObj = 
   %   cubMeanParam with properties:
   % 
   %     nInit: 300
   %      nMax: 16777216
   %
   %   
   % Author:  Fred J. Hickernell
   
   properties
      inflate %inflation factor for bounding the error
      inflateFun   % inflateFun factor 
      nInit %initial sample size
      nMax %maximum sample size
      nMu %number of integrals for solution function
      trueMuCV %true integral for control variates
   end
   
   properties (Dependent)
      nCV %number of control variates
   end
      
   properties (Hidden, SetAccess = private)
      def_nInit = 1024 %default initial number of samples
      def_nMax = 2^24 %default maximum sample size
      def_nMu = 1 %default number of integrals
      def_inflate = 1.2
      def_inflateFun = @(m) (16/3)*2.^(-m)
      def_trueMuCV = [] %default true integrals for control variates
   end
   
   
   methods
      
      % Creating a fParam process
      function obj = cubMeanParam(varargin)
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
         useDefaults = true; %true unless copying an cubMeanParam object, then false
         structInp = 0; %where is the structure
         objInp = 0; %where is the an object in the class
         if nargin %there are inputs to parse and assign
            if isa(varargin{start},'gail.cubMeanParam') 
               %the first input is a fParam object so copy it
               objInp = start;
               start = start + 1;
            end
            if nargin >= start
               if isstruct(varargin{start}) %next input is a structure containing Y
                  structInp = start;
                  start = start + 1;
               end
            end
         end
         
         done = false; %not finished parsing
         if objInp
            val = varargin{objInp}; %first input
            obj.nInit = val.nInit; %copy initial sample size
            obj.nMax = val.nMax; %copy maximum sample size
            obj.inflate = val.inflate; %copy inflation factor
            obj.nMu = val.nMu; %copy the number of means/integrals
            obj.trueMuCV = val.trueMuCV; %copy true means of control variates
            obj.inflateFun=val.inflateFun;
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
         
         f_addParamVal(p,'nInit',obj.def_nInit);
         f_addParamVal(p,'nMax',obj.def_nMax);
         f_addParamVal(p,'inflate',obj.def_inflate);
         f_addParamVal(p,'nMu',obj.def_nMu);
         f_addParamVal(p,'trueMuCV',obj.def_trueMuCV);
         f_addParamVal(p,'inflateFun',obj.def_inflateFun);
         
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
         if isfield(struct_val,'nInit')
            obj.nInit = struct_val.nInit;
         end
         if isfield(struct_val,'nMax')
            obj.nMax = struct_val.nMax;
         end
         if isfield(struct_val,'inflate')
            obj.inflate = struct_val.inflate;
         end
         if isfield(struct_val,'nMu')         
            obj.nMu = struct_val.nMu;
         end         
         if isfield(struct_val,'trueMuCV')         
            obj.trueMuCV = struct_val.trueMuCV;
         end
         
         if isfield(struct_val,'inflateFun')
            obj.inflateFun = struct_val.inflateFun;
         end
         
      end %of constructor
                            
      function set.nInit(obj,val)
         validateattributes(val, {'numeric'}, {'scalar','positive','integer'})
         obj.nInit = val;
      end
                             
      function set.nMax(obj,val)
         validateattributes(val, {'numeric'}, {'scalar','positive','integer'})
         obj.nMax = val;
      end
                                                                                      
       function set.inflate(obj,val)
         validateattributes(val, {'function_handle','numeric'}, {})
         obj.inflate = val;
       end
      
       function set.inflateFun(obj,val)
         validateattributes(val, {'function_handle','numeric'}, {})
         obj.inflateFun = val;
      end
                             
      function set.nMu(obj,val)
         validateattributes(val, {'numeric'}, {'integer', 'positive'})
         obj.nMu = val;
      end
      
      function set.trueMuCV(obj,val)
         validateattributes(val, {'numeric'}, {})
         obj.trueMuCV = val;
      end
      
      function val = get.nCV(obj)
         val = numel(obj.trueMuCV);
      end
   end
   
  methods (Access = protected)
     function propgrp = getPropertyGroups(obj)
        if ~isscalar(obj)
           propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
        else
           propList = getPropertyList(obj);
           propgrp = matlab.mixin.util.PropertyGroup(propList);
        end
     end
      
     function propList = getPropertyList(obj)
        propList = struct('nInit',obj.nInit, ...
           'nMax', obj.nMax);
        if obj.nMu ~= obj.def_nMu
           propList.nMu = obj.nMu;
        end
        if numel(obj.trueMuCV)
           propList.trueMuCV = obj.trueMuCV;
        end
     end
     
  end
   
end

