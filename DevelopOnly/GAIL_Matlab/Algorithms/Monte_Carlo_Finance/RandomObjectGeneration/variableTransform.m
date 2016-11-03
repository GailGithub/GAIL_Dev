classdef variableTransform < whiteNoise

%% variableTransform
% is a class of uniform points on different shapes.
% 

%% Properties
% This process inherits properties from the |whiteNoise| class.  Below are 
% values assigned to that are abstractly defined in that class plus some
% properties particulary for this class

   properties (SetAccess=protected) %so they can only be set by the constructor
      vtParam = struct('shape','ball') %kind of shape
   end

   properties (Constant, Hidden) %do not change & not seen
      allowShape = {'ball'} 
         %kinds of distributions that we can generate
   end

%% Methods
% The constructor for |variableTransform.| uses the |whiteNoise| constructor
% and then parses the other properties before constructing the points.

   methods
        
      % Creating a uniform shape process
      function obj = variableTransform(varargin)
         obj@whiteNoise(varargin{:}) %parse basic input       
         if isfield(obj.restInput,'vtParam')
            val = obj.restInput.vtParam;
            obj.vtParam = val;
            obj.restInput = rmfield(obj.restInput,'vtParam');
         end
      end
      
      % Set the properties of the white noise process
      function set.vtParam(obj,val)
         if isfield(val,'shape') %data for timeVector
            assert(any(strcmp(val.shape,obj.allowShape)))
            obj.vtParam.shape=val.shape;
         end
      end
                  
      % Generate uniform points on shapes
      function points=genVTPoints(obj,val)
         if strcmp(obj.vtParam.shape,'ball') && obj.timeDim.nSteps==2 %circle
            assert(strcmp(obj.wnParam.distribName,'Uniform'))
            pts = genPaths(obj,val); %get IID on square
            twopitheta = (2*pi)*pts(:,2);
            points = bsxfun(@times,sqrt(pts(:,1)), ...
               [cos(twopitheta) sin(twopitheta)]); %fix r
         end
      end
               
   end
   
   methods (Access = protected)

      function propList = getPropertyList(obj)
         propList = getPropertyList@whiteNoise(obj);
         propList.vtParam_shape = obj.vtParam.shape;
      end
      
   end

    
end

