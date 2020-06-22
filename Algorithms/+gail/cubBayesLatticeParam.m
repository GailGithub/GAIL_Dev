classdef cubBayesLatticeParam < gail.cubLatticeParam
   %GAIL.CUBBAYESLATTICEPARAM is a class containing the parameters related to
   %algorithms that find the mean of a random variable
   %   This class contains the number of integrands with the same integral,
   %   etc.
   %
   % Example 1. Construct a cubParam object with default parameters
   % >> cubBayesLatticeParamObj = gail.cubBayesLatticeParam
   % cubBayesLatticeParamObj = ***
   %
   %              f: @(x)sum(x.^2,2)
   %         domain: [2***1 double]
   %    measureType: 'uniform'
   %        measure: 'uniform'
   %         absTol: 0.0100
   %         relTol: 0
   %
   %
   % Example 2. Construct a cubParam object with properly ordered inputs
   % >> cubBayesLatticeParamObj = gail.cubBayesLatticeParam(@(x) sum(x.^3.2),[0 0; 2 2])
   % cubBayesLatticeParamObj = ***
   %
   %              f: @(x)sum(x.^3.2)
   %         domain: [2***2 double]
   %    measureType: 'uniform'
   %        measure: 'uniform'
   %         absTol: 0.0100
   %         relTol: 0
   %
   %
   % Example 3. Using name/value pairs
   % >> cubBayesLatticeParamObj = gail.cubBayesLatticeParam('domain', [0 0; 1 1], 'f', @(x) sum(x.^3.2), 'relTol', 0.1,'kerName','Ber4')
   % cubBayesLatticeParamObj = ***
   %
   %              f: @(x)sum(x.^3.2)
   %         domain: [2***2 double]
   %    measureType: 'uniform'
   %        measure: 'uniform'
   %         absTol: 0.0100
   %         relTol: 0.1000
   %        kerName: 'Ber4'
   %
   %
   % Example 4. Using a structure for input
   % >> inpStruct.f = @(x) sin(sum(x,2));
   % >> inpStruct.domain = [zeros(1,4); ones(1,4)];
   % >> inpStruct.kerName = 'Ber4';
   % >> cubBayesLatticeParamObj = gail.cubBayesLatticeParam(inpStruct)
   % cubBayesLatticeParamObj = ***
   %
   %              f: @(x)sin(sum(x,2))
   %         domain: [2***4 double]
   %    measureType: 'uniform'
   %        measure: 'uniform'
   %         absTol: 0.0100
   %         relTol: 0
   %        kerName: 'Ber4'
   %
   %
   % Example 5. Copying a cubParam object and changing some properties
   % >> NewCubBayesLatticeParamObj = gail.cubBayesLatticeParam(cubBayesLatticeParamObj,'GPMean',0)
   % NewCubBayesLatticeParamObj = ***
   %
   %              f: @(x)sin(sum(x,2))
   %         domain: [2***4 double]
   %    measureType: 'uniform'
   %        measure: 'uniform'
   %         absTol: 0.0100
   %         relTol: 0
   %        kerName: 'Ber4'
   %         GPMean: 0
   %
   %
   % Author: Fred J. Hickernell

   properties
      kerName %name of kernel function used for Bayesian cubature
      GPMean %value of the Gaussian process mean
   end

    properties (Dependent = true)
       kernel %kernel function used for Bayesian cubature
    end

   properties (Hidden, SetAccess = private)
      def_kerName = 'Ber2' %periodizing transformation
      def_GPMean = [];
   end


   methods

      % Creating a cubBayesLatticeParam process
      function obj = cubBayesLatticeParam(varargin)
         %this constructor essentially parses inputs
         %the parser will look for the following in order
         %  # a copy of a cubParam object
         %  # a structure
         %  # a function
         %  # a domain
         %  # a measure
         %  # numbers: absTol, relTol, trueMuCV, nMu, nf, inflate
         %  # name-value pairs

         start = 1; %index to begin to parse
         useDefaults = true; %true unless copying an fParam object, then false
         objInp = 0; %where is the an object in the class
         structInp = 0; %where is the structure
         if nargin %there are inputs to parse and assign
            if isa(varargin{start},'gail.cubBayesLatticeParam')
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
         end

         %Parse errorParam properties
         whichParse = [objInp structInp start:nargin];
         whichParse = whichParse(whichParse > 0);
         obj@gail.cubLatticeParam(varargin{whichParse});

         if objInp
            val = varargin{objInp}; %first input
            obj.kerName = val.kerName; %copy kernel name
            obj.GPMean = val.GPMean; %Gaussian process mean
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
         f_addParamVal(p,'kerName',obj.def_kerName);
         f_addParamVal(p,'GPMean',obj.def_GPMean);

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
         if isfield(struct_val,'kerName')
            obj.kerName = struct_val.kerName;
         end
         if isfield(struct_val,'GPMean')
            obj.GPMean = struct_val.GPMean;
         end

      end %of constructor

      function set.kerName(obj,val)
         validateattributes(val, {'char'}, {})
         obj.kerName = val;
      end

      function set.GPMean(obj,val)
         validateattributes(val, {'numeric'}, {})
         obj.GPMean = val;
      end

      function val = get.kernel(obj) %volume of the domain
         if any(strcmp(obj.kerName,{'Ber4'}))
            val = @(x,theta) prod(1 - theta*((x.*(1-x)).^2 - 1/30),2);
         else %assume kerName = Ber2
            val = @(x,theta) prod(1 + theta*(-x.*(1-x) + 1/6),2);
         end
      end

   end

   methods (Access = protected)

      function propList = getPropertyList(obj)
         propList = getPropertyList@gail.cubLatticeParam(obj);
         if ~strcmp(obj.kerName, obj.def_kerName)
            propList.kerName = obj.kerName;
         end
         if numel(obj.GPMean) ~= numel(obj.def_GPMean) || ...
            any(obj.GPMean ~= obj.def_GPMean)
            propList.GPMean = obj.GPMean;
         end
      end

   end



end
