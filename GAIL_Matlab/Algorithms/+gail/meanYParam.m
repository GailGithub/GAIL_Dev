classdef meanYParam < handle & matlab.mixin.CustomDisplay
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
   %       absTol: 0.0100
   %       relTol: 0
   %        alpha: 0.0100
   %
   %
   % Example 2. Construct a cubParam object with properly ordered inputs
   % >> meanYParamObj = gail.meanYParam(@(n) sum(rand(n,4).^3,2),0.001)
   % meanYParamObj = 
   %   meanYParam with properties:
   % 
   %            Y: @(n)sum(rand(n,4).^3,2)
   %       absTol: 1.0000e-03
   %       relTol: 0
   %        alpha: 0.0100
   %
   %
   % Example 3. Using name/value pairs
   % >> meanYParamObj = gail.meanYParam('nSig', 1e4, 'Y', @(n)sin(sum(rand(n,4).^3,2)), 'relTol', 0.1)
   % meanYParamObj = 
   %   meanYParam with properties:
   % 
   %            Y: @(n)sin(sum(rand(n,4).^3,2))
   %       absTol: 0.0100
   %       relTol: 0.1000
   %        alpha: 0.0100
   %         nSig: 10000
   %
   %
   % Example 4. Using a structure for input
   % >> inpStruct.Y = @(n) sin(sum(rand(n,2),2));
   % >> inpStruct.nSig = 1e4;
   % >> inpStruct.relTol = 0.1;
   % >> meanYParamObj = gail.meanYParam(inpStruct)
   % meanYParamObj = 
   %   meanYParam with properties:
   % 
   %            Y: @(n)sin(sum(rand(n,2),2))
   %       absTol: 0.0100
   %       relTol: 0.1000
   %        alpha: 0.0100
   %         nSig: 10000
   %
   %
   % Example 5. Copying a meanYParam object and changing some properties
   % >> NewMeanYParamObj = gail.meanYParam(meanYParamObj,'Y',@(n) rand(n,3))
   % NewMeanYParamObj = 
   %   meanYParam with properties:
   % 
   %            Y: @(n)rand(n,3)
   %       absTol: 0.0100
   %       relTol: 0.1000
   %        alpha: 0.0100
   %         nSig: 10000
   %
   %
   % Author: Fred J. Hickernell

   properties
      Y %random number generator
      alpha %uncertainty
      err %an errorParam object
      CM %a cubMeanParam object
      nSig %sample size to estimate variance
      nY %number of Y for each mean
   end
   
    properties (Dependent = true)
       nYOut %number of Y outputs
    end
   
   properties (Hidden, SetAccess = private)
      def_Y = @(n) rand(n,1) %default random number generator
      def_alpha = 0.01 %default uncertainty
      def_nSig = 1000 %default uncertainty
      def_nY = 1 %default number of Y per mean
   end
   
   
   methods
      
      % Creating a meanYParam process
      function obj = meanYParam(varargin)
         %this constructor essentially parses inputs
         %the parser will look for the following in order
         %  # a copy of a meanYParam object
         %  # a function
         %  # a structure
         %  # numbers: absTol, relTol
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
         
         if objInp
            val = varargin{objInp}; %first input
            obj.err = gail.errorParam(val.err); %copy errorParam object
            obj.CM = gail.cubMeanParam(val.CM); %copy cubMeanParam object
            obj.Y = val.Y; %copy random number generator
            obj.alpha = val.alpha; %copy uncertainty
            obj.nSig = val.nSig; %copy sample size for sigma
            obj.nY = val.nY; %copy number of Y values per mu
            useDefaults = false;
         end

         %Parse errorParam properties
         whichErrParse = [structInp start:nargin];
         whichErrParse = whichErrParse(whichErrParse > 0);
         if objInp
            obj.err = gail.errorParam(obj.err,varargin{whichErrParse});
         else
            obj.err = gail.errorParam(varargin{whichErrParse});
         end

         %Parse cubMeanParam properties
         whichCMParse = [structInp start:nargin];
         whichCMParse = whichCMParse(whichCMParse > 0);
         if objInp
            obj.CM = gail.cubMeanParam(obj.CM,varargin{whichCMParse});
         else
            obj.CM = gail.cubMeanParam(varargin{whichCMParse});
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
           parseRange = []; %to account for the two tolerances already parsed
         end
         f_addParamVal(p,'Y',obj.def_Y);
         f_addParamVal(p,'alpha',obj.def_alpha);
         f_addParamVal(p,'nSig',obj.def_nSig);
         f_addParamVal(p,'nY',obj.def_nY);
         
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
                                                    
      function set.nY(obj,val)
         validateattributes(val, {'numeric'}, {'positive','integer'})
         obj.nY = val;
      end
                                              
      function val = get.nYOut(obj)
         val = numel(obj.Y(1)); 
      end         
     
      
   end
   
   methods (Access = protected)
   
      function outval = checkNY(obj,inval)
         assert(obj.fun.nYOut - sum(inval) == obj.CM.nCV)
         outval = inval;
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
         propList = struct('Y', obj.Y, ...
            'absTol', obj.err.absTol, ...
            'relTol', obj.err.relTol, ...
            'alpha', obj.alpha);
         if obj.nSig ~= obj.def_nSig
            propList.nSig = obj.nSig;
         end
         if obj.CM.nMax ~= obj.CM.def_nMax
            propList.nMax = obj.CM.nMax;
         end
         if obj.CM.nMu ~= obj.CM.def_nMu
            propList.nMu = obj.CM.nMu;
         end
         if obj.nY ~= obj.def_nY
            propList.nY = obj.nY;
         end
         if numel(obj.CM.trueMuCV)
            propList.trueMuCV = obj.CM.trueMuCV;
         end

      end
      
   end
 
   
end

