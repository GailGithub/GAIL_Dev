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

   properties (SetAccess=public) %so they can only be set by the constructor
      %added in whiteNoise
      bmParam = struct('assembleType', 'diff') %method for assembling browniam Motion from white noise
   end

   properties (Constant, Hidden) %do not change & not seen
      allowassembleType = {'diff','PCA','BBridge'} 
   end
   
   
%% Methods
% The constructor for |brownianMotion| uses the |whiteNoise| constructor
% and then parses the other properties. The function |genBMPaths| generates
% the Brownian motion paths based on |whiteNoise| paths.

   methods
        
      % Creating a Brownian Motion process
      function obj = brownianMotion(varargin)         
         obj@whiteNoise(varargin{:}) %parse basic input
         if nargin>0
            val=varargin{1};
            if isa(val,'brownianMotion')
               obj.bmParam = val.bmParam;
               if nargin == 1
                  return
               end
            end
            if isfield(obj.restInput,'bmParam')
               val = obj.restInput.bmParam;
               obj.bmParam = val;
               obj.restInput = rmfield(obj.restInput,'bmParam');
            end
         end
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
                  sqrt(obj.timeDim.timeIncrement), ...
                  paths(:,colRange)),2);
            end
         elseif strcmp(obj.bmParam.assembleType,'PCA')
             Sigma=bsxfun(@min,obj.timeDim.timeVector',obj.timeDim.timeVector);
             [Eigenvectors,Eigenvalues]=eig(Sigma,'vector');
             [~, order] = sort(Eigenvalues, 'descend');
             A = Eigenvectors(:,order)*diag(Eigenvalues(order).^(1/2));
             paths=paths*A';
         elseif strcmp(obj.bmParam.assembleType,'BBridge')
            sobstr = sobolset(1);
            seq = sobstr(1:obj.timeDim.nSteps,1);
            [~, I] = sort(seq);
            for idx=1:obj.timeDim.dim
               colRange = ...
                  ((idx-1)*obj.timeDim.nSteps+1):idx*obj.timeDim.nSteps;
               paths(:,colRange) = cumsum(bsxfun(@times, ...
                  sqrt(obj.timeDim.timeIncrement), ...
                  paths(:,colRange(I))),2);
            end
         end
       end
                 
      function set.bmParam(obj,val)
         if isfield(val,'assembleType') %data for assemble Type
            assert(any(strcmp(val.assembleType,obj.allowassembleType)))
            obj.bmParam.assembleType = val.assembleType;
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

