classdef errorParam < handle
   %GAIL.ERRORPARAM is a class containing the parameters related to the
   %error tolerance
   %   This class contains the error tolerances, solution function, and
   %   related parameters determinign the error criterion
   %
   % Example 1.  Default values
   % >> errParamObj = gail.errorParam
   % errParamObj = 
   %   errorParam with properties:
   % 
   %       absTol: 0.010000000000000
   %       relTol: 0
   %       solFun: @(mu)mu
   %     solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
   %
   %
   % Example 2.  Assign a different absolute error tolerance
   % >> errParamObj = gail.errorParam(0.001)
   % errParamObj = 
   %   errorParam with properties:
   % 
   %       absTol: 1.000000000000000e-03
   %       relTol: 0
   %       solFun: @(mu)mu
   %     solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
   %
   %
   % Example 3. Use a name/value pair
   % >> errParamObj = gail.errorParam('relTol',0.1)
   % errParamObj = 
   %   errorParam with properties:
   % 
   %       absTol: 0.010000000000000
   %       relTol: 0.100000000000000
   %       solFun: @(mu)mu
   %     solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
   %
   %
   % Author: Fred J. Hickernell

   
   
   properties
      absTol %absolute error tolerance
      relTol %relative error tolerance
      solFun %function of vector of mu that you want to estimate
      solBdFun
      %lower and upper bounds on solFun for mu in [muhat - errbd, muhat + errbd]
   end
   
   properties (Hidden, SetAccess = private)
      def_absTol = 0.01 %default absolute error tolerance
      def_relTol = 0 %default relative error tolerance
      def_solFun = @(mu) mu %default function of vector of mu that you want to estimate
      def_solBdFun = @(muhat,errbd) [muhat - errbd, muhat + errbd]
         %lower and upper bounds on solFun for mu in [muhat - errbd, muhat + errbd]      
   end
   
   methods
      
      % Creating an errorParam object
      function obj = errorParam(varargin)
         %this constructor essentially parses inputs
         %the parser will look for the following in order
         %  # a copy of a errorParam object
         %  # a structure
         %  # numbers: absTol, relTol
         %  # name-value pairs
         
         start = 1;
         if nargin %there are inputs to parse and assign
            val=varargin{start}; %first input
            if isa(val,'gail.errorParam') %make a new copy an existing errorParam object
               obj.absTol = val.absTol; %copy absolute tolerance
               obj.relTol = val.relTol; %copy relatve tolerance
               obj.solFun = val.solFun; %copy solution function
               obj.solBdFun = val.solBdFun; %copy solution bound function
               return
            end
         end
         
         %Now begin to parse inputs
         p = inputParser; %construct an inputParser object
         p.KeepUnmatched = true; %ignore those that do not match
         p.PartialMatching = false; %don'try a partial match
         p.StructExpand = true; %expand structures
         structInp = 0; %no structure input
         done = false; %not finished parsing
         if nargin >= start
           if isstruct(varargin{start})
              structInp = start;
              start = start + 1;
           end
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
         end
         if ~done
           f_addParamVal = @addOptional;
         end
         f_addParamVal(p,'absTol',obj.def_absTol);
         f_addParamVal(p,'relTol',obj.def_relTol);
         f_addParamVal(p,'solFun',obj.def_solFun);
         f_addParamVal(p,'solBdFun',obj.def_solBdFun);
         if structInp
            parse(p,varargin{start:end},varargin{structInp}) 
            %parse inputs with a structure
         else
            parse(p,varargin{start:end}) %parse inputs w/o structure
         end
         struct_val = p.Results; %store parse inputs as a structure
         
         %Assign structure values to properties of errorParam
         obj.absTol = struct_val.absTol;
         obj.relTol = struct_val.relTol;
         obj.solFun = struct_val.solFun;
         obj.solBdFun = struct_val.solBdFun;            
      end %of constructor

      function set.absTol(obj,val)
         validateattributes(val, {'numeric'}, {'scalar','nonnegative'})
         obj.absTol = val;
      end
      
      function set.relTol(obj,val)
         validateattributes(val, {'numeric'}, {'scalar','nonnegative', ...
            '<', 1})
         obj.relTol = val;
      end
            
      function set.solFun(obj,val)
         validateattributes(val, {'function_handle'}, {})
         obj.solFun = val;
      end
            
      function set.solBdFun(obj,val)
         validateattributes(val, {'function_handle'}, {})
         obj.solBdFun = val;
      end
            
   end
   
end

