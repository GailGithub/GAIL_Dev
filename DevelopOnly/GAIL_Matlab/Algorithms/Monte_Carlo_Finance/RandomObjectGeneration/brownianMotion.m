classdef brownianMotion < whiteNoise

%% brownianMotion
% is a class of discretized stochastic processes. The marginal
% distributions are expected to be Gaussian with zero mean and variance t.
% The values of the process at different times are meant to have covariance
% equal to the minimum of the two times.
% 
%
% Example 1
% >> obj = brownianMotion
% obj = 
%   brownianMotion with properties:
% 
%                inputType: 'n'
%       timeDim_timeVector: [1 2 3]
%        timeDim_startTime: 1
%          timeDim_endTime: 3
%         timeDim_initTime: 0
%        timeDim_initValue: 0
%       wnParam_sampleKind: 'IID'
%      wnParam_distribName: 'Gaussian'
%         wnParam_xDistrib: 'Uniform'
%     bmParam_assembleType: 'diff'
%
%
% Example 2
% plot(brownianMotion)

%% Properties
% This process inherits properties from the |stochProcess| class.  Below are 
% values assigned to that are abstractly defined in that class plus some
% properties particulary for this class

   properties (SetAccess=protected) %so they can only be set by the constructor
      %added in whiteNoise
      bmParam = struct('assembleType', 'diff') %method for assembling browniam Motion from white noise
   end

%% Methods
% The constructor for |brownianMotion| uses the |whiteNoise| constructor
% and then parses the other properties. The function |genBMPaths| generates
% the Brownian motion paths based on |whiteNoise| paths.

   methods
        
      % Creating a Brownian Motion process
      function obj = brownianMotion(varargin)         
         obj@whiteNoise(varargin{:}) %parse basic input
         obj.wnParam = struct('distribName','Gaussian');
            %must have Gaussian whiteNoise paths input   
         assert(obj.timeDim.startTime >= 0) %Brownian motion goes forward in time
         obj.timeDim = struct('initTime',0,'initValue',0); 
            %Brownian motion starts at time 0 with value 0
      end
           
      % Generate Brownian Motion paths
      function paths=genPaths(obj,val)
         paths = genPaths@whiteNoise(obj,val);
         if strcmp(obj.bmParam.assembleType,'diff')
            for idx=1:obj.timeDim.dim
               colRange = ...
                  ((idx-1)*obj.timeDim.nSteps+1):idx*obj.timeDim.nSteps;
               paths(:,colRange) = cumsum(bsxfun(@times, ...
                  sqrt([obj.timeDim.timeVector(1) obj.timeDim.timeIncrement]), ...
                  paths(:,colRange)),2);
            end
         end
       end
                 
   end
   
   methods (Access = protected)

      function propList = getPropertyList(obj)
         propList = getPropertyList@whiteNoise(obj);
         propList.bmParam_assembleType = obj.bmParam.assembleType;
      end
      
   end

    
end

