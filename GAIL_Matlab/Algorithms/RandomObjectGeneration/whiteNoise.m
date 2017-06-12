classdef whiteNoise < stochProcess

%% whiteNoise
% is a class of discretized stochastic processes.
% The values of the process at different times are meant to resemble 
% independent and identically distributed (IID) random variables,
% 
% The property |paths| is a matrix corresponding to |nPaths| rows, each 
% row corresponding to a path of lenth |nSteps|.
%
% Example 1
% >> obj = whiteNoise
% obj = 
%   whiteNoise with properties:
% 
%               inputType: 'n'
%      timeDim_timeVector: [1 2 3]
%       timeDim_startTime: 1
%         timeDim_endTime: 3
%      wnParam_sampleKind: 'IID'
%     wnParam_distribName: 'Uniform'
%        wnParam_xDistrib: 'Uniform'
%
%


%% Properties
% This process inherits properties from the |stochProcess| class.  Below are 
% values assigned to that are abstractly defined in that class plus some
% properties particulary for this class

   properties (SetAccess=public)
      wnParam = struct('sampleKind','IID', ... %kind of sampling
         'distribName', 'Uniform', ... %distribution of the marginals 
         'xDistrib', 'Uniform') %kind of sampling before transformation
   end

   properties (Constant, Hidden) %do not change & not seen
      allowDistribName = {'Uniform','Gaussian'} 
         %kinds of distributions that we can generate
      allowSampleKind = {'IID','Sobol','Lattice'} 
         %kinds of sampling that we allow
      allowQRand = {'Sobol','Lattice'} 
         %kinds of samplking that we allow
   end

   properties (SetAccess=private, Hidden) %so they can only be set by the constructor
      qrandState %state of the quasi-random generator
%       defaultNPaths = 10;
%       defaultLineSpecs = {'linewidth',3}
%       defaultPointSpecs = {'markersize',25}
%       defaultFontSize = 20;
   end

%% Methods
% The constructor for |whiteNoise| uses the |stochProcess| constructor
% and then parses the other properties before constructing the paths using
% random number or quasi-random number generators.

   methods
        
      % Creating a white noise process
      function obj = whiteNoise(varargin)
         obj@stochProcess(varargin{:}) %parse basic input         
         if nargin>0
            val=varargin{1};
            if isa(val,'whiteNoise')
               obj.wnParam = val.wnParam;
               obj.qrandState = val.qrandState;
               if nargin == 1
                  return
               end
            end
            if isfield(obj.restInput,'wnParam')
               val = obj.restInput.wnParam;
               obj.wnParam = val;
               obj.restInput = rmfield(obj.restInput,'wnParam');
            end
         end
         if any(strcmp(obj.wnParam.sampleKind,obj.allowQRand)) %quasi-random numbers used
            if isfield(obj.restInput,'qrandState') 
               val = obj.restInput.qrandState; %qrandState provided
               obj.restInput = rmfield(obj.restInput,'qrandState');
            else
               val = []; %qrandState not provided
            end
            obj.qrandState = val; %set or initialize qrandstate
         end
         if strcmp(obj.inputType,'n') && ...
            strcmp(obj.wnParam.sampleKind,'IID') && ... %easier to sample from randn
            strcmp(obj.wnParam.distribName,'Gaussian') %if you want Gaussian
            obj.wnParam.xDistrib = 'Gaussian';
         end
      end
      
      % Set the properties of the white noise process
      function set.wnParam(obj,val)
         if isfield(val,'sampleKind') %data for timeVector
            assert(any(strcmp(val.sampleKind,obj.allowSampleKind)))
            obj.wnParam.sampleKind=val.sampleKind; %row
         end
         if isfield(val,'distribName') %distribName is provided
            assert(any(strcmp(val.distribName,obj.allowDistribName)))
            obj.wnParam.distribName = val.distribName;
         end
         if isfield(val,'xDistrib') %xDistrib is provided
            assert(any(strcmp(val.xDistrib,obj.allowDistribName)))
            obj.wnParam.xDistrib = val.xDistrib;
         elseif ~isfield(obj.wnParam,'xDistrib')
            obj.wnParam.xDistrib = obj.allowDistribName{1};
         end
         if isfield(val,'sampleKind') %sampleKind is provided
            assert(any(strcmp(val.sampleKind,obj.allowSampleKind)))
%             assert(strcmp(obj.wnParam.xDistrib,'Uniform') || ...
%                ~any(strcmp(obj.sampleKind,obj.allowQRand)))
            obj.wnParam.sampleKind = val.sampleKind;
         end
      end
      
      function set.qrandState(obj,val)
         if ~isempty(val) %state of qrand provided
            assert(size(val.PointSet,2) == obj.timeDim.nSteps);
               %size of point set must match nSteps
            obj.qrandState = val;
         else
            obj.qrandState = qrandstream(scramble( ...
               sobolset(obj.timeDim.nSteps),'MatousekAffineOwen'));
         end
      end
            
      % Generate white noise paths
      function paths=genPaths(obj,val)
         if strcmp(obj.inputType,'n')
            %construct paths from random number generators
            validateattributes(val,{'numeric'},{'scalar','positive'})
            if strcmp(obj.wnParam.sampleKind,'IID') %IID samples
               if strcmp(obj.wnParam.xDistrib,'Uniform') %uniform IID
                  paths=rand(val,obj.timeDim.nCols);
               elseif strcmp(obj.wnParam.xDistrib,'Gaussian') %Gaussian IID
                  paths=randn(val,obj.timeDim.nCols);
               end
            elseif strcmp(obj.wnParam.sampleKind,'Sobol') %Sobol samples
               paths=rand(obj.qrandState,val,obj.timeDim.nCols); 
                  %uniform Sobol
            end
         else
            validateattributes(val,{'numeric'},{'2d'})
            assert(size(val,2)==obj.timeDim.nCols, ...
               ['# of columns of ''x'' input to genPaths must equal' ...
               ' (# of times)(dimension)'])
            paths=val;
         end
         
         if strcmp(obj.wnParam.xDistrib,'Uniform') && ...
               strcmp(obj.wnParam.distribName,'Gaussian') %need a transformation
            paths=gail.stdnorminv(paths);
         end
      end
      
%       function varargout = plot(obj,varargin)
%          assert(strcmp(obj.inputType,'n'), ...
%             'plot requires inputType to be ''n''')
%          offset = 0;
%          if numel(varargin) >= 1
%             if any(strcmp(varargin{1},{'line','point'}))
%                plotKind = varargin{1};
%                offset = 1;
%             else
%                plotKind = 'line';
%             end
%          else
%             plotKind = 'line';
%          end
%          if numel(varargin) > offset
%             nPaths = varargin{offset+1};
%          else
%             nPaths = obj.defaultNPaths; %default 
%          end
%          paths = genPaths(obj,nPaths);
%          if strcmp(plotKind,'point')
%             if obj.timeDim.nSteps >= 2;
%                h = plot(paths(:,1),paths(:,2),'.');
%                if numel(varargin) > offset + 1
%                   set(h,varargin{offset+2:end});
%                else
%                   set(h,obj.defaultPointSpecs{:});
%                end
%             else
%                plotKind = 'line';
%             end
%          end
%          if strcmp(plotKind,'line')
%             timeVec = obj.timeDim.timeVector;
%             if numel(obj.timeDim.initTime)
%                timeVec = [obj.timeDim.initTime timeVec];
%                paths = [repmat(obj.timeDim.initValue,nPaths,1) paths];
%             end
%             h = plot(timeVec,paths,'-');
%             if numel(varargin) > offset + 1
%                set(h,varargin{offset+2:end});
%             else
%                set(h,obj.defaultLineSpecs{:});
%             end
%          end
%          set(gca,'fontsize',obj.defaultFontSize)
%          if nargout
%             varargout{1}=h;
%          end
%       end
          
   end
   
   methods (Access = protected)

      function propList = getPropertyList(obj)
         propList = getPropertyList@stochProcess(obj);
         propList.wnParam_sampleKind = obj.wnParam.sampleKind;
         propList.wnParam_distribName = obj.wnParam.distribName;
         propList.wnParam_xDistrib = obj.wnParam.xDistrib;
      end
      
   end

    
end

