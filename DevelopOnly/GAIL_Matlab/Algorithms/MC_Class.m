classdef MC_Class < handle
    %MC_Class for testing arrays
    
    properties
        randGen %rand number generator
        n %number of samples
    end
    
    methods
        function obj = MC_Class(varargin)
            obj.randGen = varargin{1};
            obj.n = varargin{2};
        end
        
        function val = genMu(obj)
            nObj = numel(obj);
            val(nObj) = 0;
            for i = 1:numel(obj)
              val(i) = mean(obj(i).randGen(obj(i).n));
            end
        end
    end
    
    
end

