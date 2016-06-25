classdef MC_Class < handle
    %MC_Class for testing arrays    
    properties
        randGen %rand number generator
        n %number of samples
    end
    
    methods
        function obj = MC_Class(randGen, n)
            if nargin ~= 0
            len = length(randGen);
            obj(1,len) = MC_Class;
            for i = 1:len
                obj(i).randGen = randGen{i};
                obj(i).n = n{i};
            end
            end
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

