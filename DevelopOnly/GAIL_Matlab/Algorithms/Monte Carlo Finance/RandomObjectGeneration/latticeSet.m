classdef latticeSet<handle
%% latticeSet
% This class let you create an  latticeSet object.
%
% Example 1
% obj = latticeSet
% obj = 
%   latticeSet with properties:
% 
%     genVector: [1x100 double]  %generating vector
%           dim: 2               %number of dimensions
%         shift: 1               %allows shift or not.
%         index: 1               %point in which sample starts.
%
% With such object you are allowed to use the method:
% 
% * genLattice(obj,N): Returns a matrix cointaing N points of the Lattice
% (Different Lines are different points. Different collumns are different
% coordinates)

    properties(SetAccess = 'protected')
       %genVector:= generating vector
        genVector = [1, 27808567, 14629093, 17951525, 873611, 20915065, 21016621, 23276977,...
        13573889, 23031381, 30852487, 26266751, 28027657, 5189743, 16396083,...
        17119083, 13327883, 5182645, 736279, 16190439, 4281341, 18011305,...
        26892147, 2860863, 18897101, 20241879, 31482589, 9705235, 23408339,...
        14422975, 27614503, 19560933, 1623155, 25844767, 4275455, 19908921,...
        32077639, 13483691, 7500601, 23724921, 20105911, 10850327, 12463705,...
        20419927, 15450237, 13699327, 5391239, 20988989, 2292675, 22914301,...
        11579857, 3631277, 29192369, 27493855, 26733513, 18979447, 7895911,...
        28510509, 5136563, 6010773, 22244255, 18318415, 21262503, 2703747,...
        12557745, 29189357, 24888579, 5272959, 18244419, 10852353, 19495673,...
        15385229, 17652639, 21737265, 13726991, 14542113, 16618901, 1921655,...
        29004257, 24806569, 2515121, 25540279, 3051655, 5637047, 19239479,...
        25413861, 19026129, 7717595, 16939323, 5596877, 3411409, 29006229,...
        29289975, 17980543, 25987481, 23186729, 4858187, 25345023, 30536647, 865955]; 
       %dim := number of dimensions:
        dim = 2;
       %shift := add a normal random vector to all points in the sample.
        shift = true; 
       %index := starting point
        index=0;
    end
    
    properties (Dependent = true)
        shiftValue 
    end
    
    
    methods
        
        function obj = latticeSet(genVector,dim,shift,index)%Constructor, with data treating
            if nargin == 0
            elseif size(genVector,1)~=1 && size(genVector,2)~=1 && length(genVector)<dim
                disp('genVector must be a vector with length smaller than dim')
            elseif length(dim)~=1 || dim<1 || length(index)~=1 || index<1
                disp('dim and index must be numbers greater than 1')
            elseif shift~=0 && shift~=1
                disp('shift must be a boolean value')
            else
              obj.genVector = genVector;
              obj.dim = dim;
              obj.shift = shift;
              obj.index = vanDerCorput(index);
            end
            obj.shiftValue = rand(1,obj.dim);
        end
        
        
        function lattice = genLattice(obj,N)%Generates the Lattice set based the latticeSet object and the amount of points required
            assert(obj.dim<=numel(obj.genVector), ...
               'dim must be no greater than the number of elements in the generating vector')
            obj.index = vanDerCorput(obj.index,length(obj.index)+N)
            lattice = obj.index(length(obj.index)-N+1:length(obj.index))*obj.genVector(1:obj.dim);
            if obj.shift
                lattice = bsxfun(@plus,lattice,obj.shiftValue);
            end
            lattice = mod(lattice,1);
        end
        
        function set.dim(obj,val)
            validateattributes(val,{'numeric'},{'scalar','positive','integer'})
            obj.dim = val;
        end 
        
        function set.shiftValue(obj,~)
            obj.shiftValue = rand(1,obj.dim);
        end 
        
        function set.genVector(obj,val)
            validateattributes(val,{'numeric'},{'vector','integer'})
            obj.genVector=val;
        end 
        
        function set.shift(obj,val)
            validateattributes(val,{'logical'},{'scalar'})
            obj.shift=val;
            
        end
            
    end
 
end

