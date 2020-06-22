classdef stochProcess < handle & matlab.mixin.CustomDisplay

%% stochProcess
% is an abstract class of discretized stochastic processes.
%
%   Concrete subclasses include 
%      o whiteNoise
%      o brownianMotion
%
% Example 1
% >> obj = stochProcess
% obj = 
%   stochProcess with properties:
% 
%              inputType: 'n'
%     timeDim_timeVector: [1 2 3]
%      timeDim_startTime: 1
%        timeDim_endTime: 3
%            timeDim_dim: 1
%
%
% Authors: Fred J. Hickernell

           
   properties
      timeDim = struct( ...
         'timeVector', 1:3, ... %vector of increasing times for process 
         'startTime', 1, ... %starting time, usually larger than zero
         'endTime', 3, ... %ending time
         'nSteps', 3, ... %number of different times
         'timeIncrement', [1 1], ... %increment between time intervals
         'dim', 1, ... %dimension of process
         'nCols', 3, ... %number of columns of the process matrix = nSteps x dim
         'initTime', [], ... %initial time, normally zero if it exists
         'initValue', []) %initial value, if it exists
      inputType = 'n' %input type: 'n' for number of paths, 
                      %            'x' for array of numbers
   end
   
   properties (Hidden, SetAccess=protected)
      restInput = []; %remaining input not yet parsed
   end

   properties (Hidden, SetAccess=private) %so they can only be set by the constructor
      defaultNPaths = 10
      defaultSpecs = {'linewidth',3,'markersize',25}
      defaultFontSize = 20
      defaultPlotKind = 'yt-'
      allowedPlotKind = {'yt.','yt-','yy','hist'}
      defaultColor = [0 0.447 0.741]; %MATLAB blue
   end
   

   methods
        
      % Creating a stochastic process
      function obj = stochProcess(varargin)
         %this constructor essentially parses inputs
 
         if nargin %there are inputs to parse and assign
            val=varargin{1}; %first input
            if isa(val,'stochProcess') %make a copy an existing stochProcess object
               obj.inputType = val.inputType;
               obj.timeDim = val.timeDim;
               if nargin == 1
                   return
               else
                   val = varargin{2}; %structure that can contain additional inputs                  
               end
            end
            obj.restInput = val; %rest of inputs for subclasses
            if isstruct(val) %structure with properties of stochProcess
               if isfield(val,'inputType')
                  obj.inputType = obj.restInput.inputType;
                  obj.restInput = rmfield(obj.restInput,'inputType');
               end
               if isfield(val,'timeDim')
                  obj.timeDim = val.timeDim;
                  obj.restInput = rmfield(obj.restInput,'timeDim');
               end
            end
         end
      end
      
      function set.timeDim(obj,val) %set the timeDim property
         if isfield(val,'timeVector') %data for timeVector provided
            obj.timeDim.timeVector = val.timeVector(:)'; %row vector
              MATLABVERSION = gail.matlab_version;
              if MATLABVERSION <= 8.1
                  if any(diff(obj.timeDim.timeVector)<=0)
                       error('Need obj.timeDim.timeVector with increasing values.')               
                  end
              else
                  validateattributes(obj.timeDim.timeVector, ...
                      {'numeric'}, {'increasing'})
              end
            obj.timeDim.nSteps = numel(obj.timeDim.timeVector); %number of time steps
         elseif isfield(val,'nSteps') %data for nSteps provided
            validateattributes(val.nSteps, ...
               {'numeric'}, {'scalar','integer','positive'})
            obj.timeDim.timeVector = 1:val.nSteps; %time vector with intervals of one    
         end
         if isfield(val,'dim')
            validateattributes(val.dim, {'numeric'}, ...
               {'scalar','integer','positive'})
            obj.timeDim.dim = val.dim; %dimension of the stochastic process            
         end
         if isfield(val,'initTime')
            if numel(val.initTime)
               validateattributes(val.initTime, {'numeric'},{'scalar'})
               assert(val.initTime <= obj.timeDim.startTime)
               obj.timeDim.initTime = val.initTime; %initial time before startTime
            end
         end
         if isfield(val,'initValue')
            if numel(val.initTime)
               validateattributes(val.initValue, {'numeric'},{'vector'})
               obj.timeDim.initValue = val.initValue; %initial value of the stochastic process
            end
         end
         %compute all of the dependent properties
         obj.timeDim.startTime = obj.timeDim.timeVector(1); %start time
         obj.timeDim.endTime = obj.timeDim.timeVector(end); %end time
         obj.timeDim.nSteps = numel(obj.timeDim.timeVector); %number of steps
         obj.timeDim.timeIncrement = diff([obj.timeDim.initTime obj.timeDim.timeVector]); %increment between time steps
         obj.timeDim.nCols = obj.timeDim.nSteps * obj.timeDim.dim; %total number of columns
      end
            
      function set.inputType(obj,val)
         if strcmp(val,'x') %paths are input
            obj.inputType = 'x'; %input matrices of paths
         else
            obj.inputType = 'n'; %input numbers of paths
         end
      end
      
      function varargout = plot(obj,varargin) %%DOCUMENT & FIX
         if ~any(strcmp(methods(obj),'genPaths'))
            warning(['Objects of class ' class(obj) ' cannot be plotted' ...
               ' because no genPaths method exists'])
            return %can only plot if a genPaths method exists
         end
         offset = 1; %the last input just looked at 
         if nargin > offset %another input exists
            if ischar(varargin{offset})
               if any(strcmp(varargin{offset},obj.allowedPlotKind)) %is it a kind of plot
                  plotKind = varargin{offset}; %then set it to the desired kind
                  offset = offset+1; %increase offset
               else
                  warning('Second input is not a valid kind of plot')
                  return
               end
            else
               plotKind = obj.defaultPlotKind; %otherwise set to default
            end
         else
            plotKind = obj.defaultPlotKind; %otherwise set plotKind to default
         end
         if nargin > offset %does another input exist
            genPathsInput = varargin{offset}; %input to the method to generate paths
            offset = offset + 1;
         elseif strcmp(obj.inputType,'n')
            genPathsInput = obj.defaultNPaths; %default for # of paths to be plotted
         else
            error('Need input for inputType ''x''')               
         end
         paths = genPaths(obj,genPathsInput); %generate paths to be plotted
         if strcmp(plotKind,'yy') %plot path columns against each other
            if obj.timeDim.nSteps >= 2
               h = plot(paths(:,1),paths(:,2),'.');
            else
               plotKind = obj.defaultPlotKind;
            end
         end
         if strncmp(plotKind,'yt',2) || strcmp(plotKind,'hist')
            timeVec = obj.timeDim.timeVector; %times to be plotted at
            nPaths = size(paths,1); %number of paths to be plotted
            if numel(obj.timeDim.initTime)
               timeVec = [obj.timeDim.initTime timeVec]; %add initial time to time vector               
               paths = [repmat(obj.timeDim.initValue,nPaths,1) paths]; %add initial value to paths
            end
            if strncmp(plotKind,'yt',2) %a y versus time plot
               h = plot(timeVec,paths,plotKind(3)); %the third argument gives line type
            else %a histogram
               nTimeVec = numel(timeVec);
               maxHtVec = [diff(timeVec) 0];
               maxHtVec(nTimeVec) = maxHtVec(nTimeVec-1);
               ptsPerBin = ceil(sqrt(nPaths));
               binBreaks = (1:ptsPerBin:nPaths)';
               h = zeros(1,nTimeVec);
               for ii = 1:nTimeVec
                  pathVals = sort(paths(:,ii));        
                  binBrkVals = pathVals(binBreaks);
                  binWdths = diff(binBrkVals);
                  binCts = diff(binBreaks);
                  binCts(1) = binCts(1)+1;
                  tempa = min(binWdths./binCts);
                  tempb = (3 + 0.2*nPaths*tempa)/(4 + 0.8*nPaths*tempa);
                  binHts = tempb*ones(size(binCts))*maxHtVec(ii);
                  bWNotZero = binWdths>0;
                  binHts(bWNotZero) = (tempa*tempb*maxHtVec(ii))* ...
                     (binCts(bWNotZero)./binWdths(bWNotZero));
                  patchy = repmat(binBrkVals',2,1);
                  patchy = [patchy(:); patchy(1)];
                  patchx = [0 binHts'; binHts' 0];
                  patchx = [patchx(:); 0];
                  h(ii) = patch(patchx+timeVec(ii),patchy,obj.defaultColor); 
                  hold on
               end
               set(h,'EdgeColor',obj.defaultColor);
            end                
         end
         if nargin > offset %adjust marker and line sizes
            set(h,varargin{offset:end});
         else
            set(h,obj.defaultSpecs{:});
         end
         set(gca,'fontsize',obj.defaultFontSize) %set font size
         if nargout
            varargout{1} = h; %output plot handle
         end
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
      propList = struct('inputType',obj.inputType, ...
         'timeDim_timeVector', obj.timeDim.timeVector, ...
         'timeDim_startTime', obj.timeDim.startTime, ...
         'timeDim_endTime', obj.timeDim.endTime);
      if numel(obj.timeDim.initTime)
         propList.timeDim_initTime = obj.timeDim.initTime;
      end
      if numel(obj.timeDim.initValue)
         propList.timeDim_initValue = obj.timeDim.initValue;
      end
      propList.timeDim_dim = obj.timeDim.dim;
   end

   end
          
end

