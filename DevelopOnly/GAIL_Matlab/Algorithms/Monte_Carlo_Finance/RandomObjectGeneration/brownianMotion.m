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
      allowassembleType = {'diff','PCA','bridge','bridgeApx'} 
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
                  sqrt([obj.timeDim.timeVector(1) obj.timeDim.timeIncrement]), ...
                  paths(:,colRange)),2);
            end
         elseif strcmp(obj.bmParam.assembleType,'PCA')
             Sigma=bsxfun(@min,obj.timeDim.timeVector',obj.timeDim.timeVector);
             [Eigenvectors,Eigenvalues]=eig(Sigma);
             A = Eigenvectors*Eigenvalues.^(1/2);
             paths=paths*A';
         elseif strcmp(obj.bmParam.assembleType,'bridge')
             t = repmat([0,obj.timeDim.timeVector],val,1); %insert t0 and replicate matrix for val samples 
             itrlmt=floor(log2(obj.timeDim.nSteps)); %set the maximum iteration 
             for idx=1:obj.timeDim.dim
                 basepaths = [zeros(val,1),...
                     paths(:,((idx-1)*obj.timeDim.nSteps+1):idx*obj.timeDim.nSteps)]; %create basepaths for each dimemsion and insert B(0)
                 basepaths(:,end) = sqrt(t(:,end)).*basepaths(:,end); %insert B(T) to the basepaths
                 insertedIdx = [1,obj.timeDim.nSteps+1]; %create an array record the index of inserted time 
                 for i = 1:itrlmt
                     currIdx = floor((insertedIdx(1:end-1)+insertedIdx(2:end))./2); %calculate the midpoint of previous time interval to be inserted
                     basepaths(:,currIdx )= (basepaths(:,insertedIdx(1:end-1)).*...
                         (t(:,insertedIdx(2:end))-t(:,currIdx))+basepaths(:,insertedIdx(2:end)).*...
                         (t(:,currIdx) - t(:,insertedIdx(1:end-1))))./(t(:,insertedIdx(2:end))-...
                         t(:,insertedIdx(1:end-1)))+sqrt((t(:,insertedIdx(2:end))-t(:,currIdx)).*...
                         (t(:,currIdx) - t(:,insertedIdx(1:end-1)))./(t(:,insertedIdx(2:end))-...
                         t(:,insertedIdx(1:end-1)))).* paths(:,currIdx-1+(idx-1)*obj.timeDim.nSteps);
                     %insert the brownian motion to the basepaths
                     insertedIdx = sort([insertedIdx, currIdx]);
                     %update insertedIdx with the index just added
                 end
                 if obj.timeDim.nSteps~=2^itrlmt %check whether all time inserted
                    missedIdx = setdiff(1:obj.timeDim.nSteps+1,insertedIdx); %find the missing time
                    basepaths(:,missedIdx) = (basepaths(:,missedIdx-1).*...
                         (t(:,missedIdx+1)-t(:,missedIdx))+basepaths(:,missedIdx+1).*...
                         (t(:,missedIdx) - t(:,missedIdx-1)))./(t(:,missedIdx+1)-...
                         t(:,missedIdx-1))+sqrt((t(:,missedIdx+1)-t(:,missedIdx)).*...
                         (t(:,missedIdx) - t(:,missedIdx-1))./(t(:,missedIdx+1)-...
                         t(:,missedIdx-1))).* paths(:,missedIdx-1+(idx-1)*obj.timeDim.nSteps);
                    %insert the rest brownian motion
                 end
                 paths(:,((idx-1)*obj.timeDim.nSteps+1):idx*obj.timeDim.nSteps) = basepaths(:,2:end);
                 % assign value of basepaths into paths according to different dimemsion
                 % skip B(0)
             end
         
         elseif strcmp(obj.bmParam.assembleType,'bridgeApx')
             X = paths; 
             mk = floor(log2(1:obj.timeDim.nSteps-1));
             p = sobolset(1);
             tkOvT = 1 - p(2:obj.timeDim.nSteps)';
             for idx=1:obj.timeDim.dim
                 colRange = ...
                  ((idx-1)*obj.timeDim.nSteps+1):idx*obj.timeDim.nSteps;
                 for t_idx=1:obj.timeDim.nSteps
                     paths(:,t_idx+(idx-1)*obj.timeDim.nSteps) = X(:,colRange(1)).*...
                         obj.timeDim.timeVector(t_idx)./sqrt(obj.timeDim.endTime) + ...
                         sum( X(:,colRange(2:end)).*repmat(sqrt(obj.timeDim.endTime./...
                         (2.^(mk+2))).*(1-min(abs((2.^(mk+1)).*(obj.timeDim.timeVector(t_idx)./...
                         obj.timeDim.endTime-tkOvT)),1)),size(paths,1),1),2);
                  end
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

