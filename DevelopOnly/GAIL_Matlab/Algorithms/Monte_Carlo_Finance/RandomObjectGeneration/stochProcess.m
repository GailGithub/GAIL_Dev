classdef stochProcess < handle & matlab.mixin.CustomDisplay
% is a an abstract class of discretized stochastic processes.
%   Concrete subclasses include 
%      o white noise
%      o Brownian motion
%
% Example 1
% >> stochProcess
% ans = 
%   stochProcess with properties:
% 
%              inputType: 'n'
%     timeDim_timeVector: [1 2 3]
%      timeDim_startTime: 1
%        timeDim_endTime: 3
%
%
% Authors: Fred J. Hickernell

           
   properties
      timeDim = struct( ...
         'timeVector', 1:3, ... %vector of times where process is discretized
         'startTime', 1, ... %starting time
         'endTime', 3, ... %ending time
         'nSteps', 3, ... %number of different times
         'timeIncrement', [1 1], ... %increment between time intervals
         'dim', 1, ... %dimension of process
         'nCols', 3, ... %number of columns of the process matrix = nSteps*dim
         'initTime', [], ... %initial time, normally zero if it exists
         'initValue', []) %initial value
      inputType = 'n' %input type: 'n' for number of paths, 
                      %            'x' for array of numbers
   end
   
   properties (Hidden, SetAccess=protected)
      restInput = []; %remaining input not yet parsed
   end

   methods
        
      % Creating a stochastic process
      function obj = stochProcess(varargin)
         %this constructor essentially parses inputs
         
         if nargin>0
            val=varargin{1};
            if isa(val,'stochProcess')
               obj.inputType = val.inputType;
               obj.timeDim = val.timeDim;
               if nargin == 1
                   return
               else
                   val = varargin{2};
               end
            end
            if isstruct(val)
               obj.restInput = val;
               if isfield(val,'inputType')
                  obj.inputType = val.inputType;
                  obj.restInput=rmfield(obj.restInput,'inputType');
               end
               if isfield(val,'timeDim')
                  obj.timeDim = val.timeDim;
                  obj.restInput=rmfield(obj.restInput,'timeDim');
               end
            end
         end
      end
      
      function set.timeDim(obj,val) %set the timeDim property
         if isfield(val,'timeVector') %data for timeVector
            obj.timeDim.timeVector = val.timeVector(:)'; %row
            validateattributes(obj.timeDim.timeVector, ...
               {'numeric'}, {'increasing'})
            obj.timeDim.nSteps = numel(obj.timeDim.timeVector); %number of steps
         elseif isfield(val,'nSteps') %data for nSteps provided
            validateattributes(val.nSteps, ...
               {'numeric'}, {'scalar','integer','positive'})
            obj.timeDim.timeVector = 1:obj.timeDim.nSteps; %time vector      
         end
         if isfield(val,'dim')
            validateattributes(val.dim, {'numeric'}, ...
               {'scalar','integer','positive'})
            obj.timeDim.dim=val.dim; %dimension
         end
         if isfield(val,'initTime')
            if numel(val.initTime)>0 
               validateattributes(val.initTime, {'numeric'},{'scalar'})
               assert(val.initTime <= obj.timeDim.startTime)
               obj.timeDim.initTime=val.initTime; %initial time before startTime
            end
         end
         if isfield(val,'initValue')
            if numel(val.initTime)>0 
               validateattributes(val.initValue, {'numeric'},{'vector'})
               obj.timeDim.initValue=val.initValue; %initial value
            end
         end
         %compute all of the dependent properties
         obj.timeDim.startTime = obj.timeDim.timeVector(1); %start time
         obj.timeDim.endTime = obj.timeDim.timeVector(end); %end time
         obj.timeDim.nSteps = numel(obj.timeDim.timeVector); %number of steps
         obj.timeDim.timeIncrement = diff(obj.timeDim.timeVector); %increment between time steps
         obj.timeDim.nCols = obj.timeDim.nSteps*obj.timeDim.dim; %total number of columns
      end
            
      function set.inputType(obj,val)
         if strcmp(val,'x') %paths are input
            obj.inputType = 'x'; %input matrices of paths
         else
            obj.inputType = 'n'; %input numbers of paths
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
   end
   

end
       
end

