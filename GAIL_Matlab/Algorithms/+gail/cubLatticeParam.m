classdef cubLatticeParam < gail.cubParam
   %GAIL.CUBLATTICEPARAM is a class containing the parameters related to
   %algorithms that find the mean of a random variable
   %   This class contains the number of integrands with the same integral,
   %   etc.
   %
   % Example 1. Construct a cubParam object with default parameters
   % >> cubLatticeParamObj = gail.cubLatticeParam
   % cubLatticeParamObj = ***
   %
   %              f: @(x)sum(x.^2,2)
   %         domain: [2***1 double]
   %    measureType: 'uniform'
   %        measure: 'uniform'
   %         absTol: 0.0100
   %         relTol: 0
   %
   %
   % Example 2. Using name/value pairs
   % >> cubLatticeParamObj = gail.cubLatticeParam('domain', [-2 -2; 2 2], 'f', @(x) sum(x.^3.2), 'relTol', 0.1, 'isShift', false)
   % cubLatticeParamObj = ***
   %
   %              f: @(x)sum(x.^3.2)
   %         domain: [2***2 double]
   %    measureType: 'uniform'
   %        measure: 'uniform'
   %         absTol: 0.0100
   %         relTol: 0.1000
   %        isShift: 0
   %
   %
   % Example 3. Using a structure for input
   % >> inpStruct.f = @(x) sin(sum(x,2));
   % >> inpStruct.domain = [zeros(1,4); ones(1,4)];
   % >> inpStruct.isShift = false;
   % >> cubLatticeParamObj = gail.cubLatticeParam(inpStruct)
   % cubLatticeParamObj = ***
   %
   %              f: @(x)sin(sum(x,2))
   %         domain: [2***4 double]
   %    measureType: 'uniform'
   %        measure: 'uniform'
   %         absTol: 0.0100
   %         relTol: 0
   %        isShift: 0
   %
   %
   % Example 4. Copying a cubParam object and changing some properties
   % >> NewCubLatticeParamObj = gail.cubLatticeParam(cubLatticeParamObj,'measure','Lebesgue')
   % NewCubLatticeParamObj = ***
   %
   %              f: @(x)sin(sum(x,2))
   %         domain: [2***4 double]
   %    measureType: 'uniform'
   %        measure: 'Lebesgue'
   %         absTol: 0.0100
   %         relTol: 0
   %        isShift: 0
   %
   %
   % Author: Fred J. Hickernell

   properties
      periodTransform %periodizing transformation
      isShift %is the lattice shifted
   end

   properties (Dependent = true)
      fff %function after periodizing transformation
   end


   properties (Hidden, SetAccess = private)
      def_periodTransform = 'tent' %default periodizing transformation
      def_isShift = true %default is a random shift
   end


   methods

      % Creating a cubParam process
      function obj = cubLatticeParam(varargin)
         %this constructor essentially parses inputs

         start = 1; %index to begin to parse
         useDefaults = true; %true unless copying an fParam object, then false
         objInp = 0; %where is the an object in the class
         structInp = 0; %where is the structure
         if nargin %there are inputs to parse and assign
            if isa(varargin{start},'gail.cubLatticeParam')
               %the first input is a cubLatticeParam object so copy it
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

         %Parse errorParam properties
         whichParse = [objInp structInp start:nargin];
         whichParse = whichParse(whichParse > 0);
         obj@gail.cubParam(varargin{whichParse});

         if objInp
            val = varargin{objInp}; %first input
            obj.periodTransform = val.periodTransform; %copy integration measure
            obj.isShift = val.isShift; %copy whether to shift
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
               parseRange = start:nargin;
               done = true;
            end
         end
         if ~done %then nothingleft or just numbers
            f_addParamVal = @addOptional;
            parseRange = []; %to account for the two tolerances already parsed
         end
         f_addParamVal(p,'periodTransform',obj.def_periodTransform);
         f_addParamVal(p,'isShift',obj.def_isShift);

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
         if isfield(struct_val,'periodTransform')
            obj.periodTransform = struct_val.periodTransform;
         end

         %Assign values of structure to corresponding class properties
         if isfield(struct_val,'isShift')
            obj.isShift = struct_val.isShift;
         end

      end %of constructor

      function set.periodTransform(obj,val)
         validateattributes(val, {'char'}, {})
         obj.periodTransform = val;
      end

      function set.isShift(obj,val)
         validateattributes(val, {'logical'}, {'scalar'})
         obj.isShift = val;
      end

      function val = get.fff(obj)
         if strcmp(obj.periodTransform,'tent') && strcmp(obj.domainType,'box')
            val = @(x) obj.ff(1 - abs(2*x-1)); %tent transformation
         else %no transformation
            val = obj.ff;
         end

      end


   end

   methods (Access = protected)

      function propList = getPropertyList(obj)
         propList = getPropertyList@gail.cubParam(obj);
         if ~strcmp(obj.periodTransform,obj.def_periodTransform)
            propList.periodTransform = obj.periodTransform;
         end
         if obj.isShift ~= obj.def_isShift
            propList.isShift = obj.isShift;
         end
      end

   end



end
