classdef latticeSet<handle
%% latticeSet
% This class let you create an  latticeSet object.
%
% Example 1
%
% obj = latticeSet
%
% obj = 
%
%   latticeSet with properties:
% 
%     genVector: [1x100 double]  %generating vector
%           dim: 2               %number of dimensions
%         shift: 1               %allows shift or not.
%    shiftValue: [0.8147 0.9058] %the vector that shifts every point
%         index: 0               %point in which sample starts.
%
% With such object you are allowed to use the method genLattice(obj,N), which 
% returns a matrix cointaing N points of the Lattice
% (Different Lines are different points. Different collumns are different
% coordinates)
%
% Authors: Maurício Darós Andrade, Lluís Antoni Jiménez Rugama.


%% Properties
%
% We start with the properties:
    properties
       %shift := add a normal random vector to all points in the sample.
        shift = true;
        shiftValue; 
        %index := starting point
        index=0;
    end

%%
% For properties that makes change on other properties when they are changed
% we need to set them as dependent and create auxiliar private properties.
% Here we have dim and genVector as dependent properties. Also we have
% privateDim and privateGenVector. These last ones will carry the value
% wanted for dim and genVector. Whenever we want to access or change dim or
% genVector we will use get and set methods. These dependent properties are
% only there to serve as label for our true properties (that is, the private ones
% that carry the values). Furthermore, the dependent properties allow us to
% change other properties when changing them, what is not possible with the
% private properties.
    properties(Dependent=true)
       %dim := number of dimensions:
        dim
        %genVector:= generating vector
        genVector
    end
    
    properties(Access=private)
       privateDim=2;
       privateGenVector = [1, 27808567, 14629093, 17951525, 873611, 20915065, 21016621, 23276977,...
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

    end
    
%% Methods 
% We have the constructor method which assign values to the
% properties upon some treating.
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
              obj.index = index;
            end
            obj.shiftValue = rand(1,obj.dim);
        end
%%
%
% The following method generates the Lattice set based the latticeSet object
% and the amount of points required. It will generate points starting from
% the current index and then change the index to index+N, so when you call
% the function more than one time you don't get repeated points.
% 
% It basically works using the function vanDerCorput, that returns us N Van
% Der Corput points. It takes these points and multiplies by the generating
% vector. After doing this, it adds to the points an optional shift.
% Finally, it takes all coordinates modulo 1.
%
        function lattice = genLattice(obj,N)
            lattice = vanDerCorput(obj.index,obj.index+N-1)*obj.genVector(1:obj.dim);
            if obj.shift
                lattice = bsxfun(@plus,lattice,obj.shiftValue);
            end
            lattice = mod(lattice,1);
            obj.index = obj.index+N;
        end
        
%% Get Methods
% These allows you to access dim and genVector by getting their values
% on the properties privateDim and privateGenVector.

        function dim = get.dim(obj)
            dim=obj.privateDim;
        end
        function genVector = get.genVector(obj)
            genVector=obj.privateGenVector;
        end
        
%% Set methods
% These set functions treat any new values assigned to the properties and
% work together with other properties so you don't have to worry about all
% properties. The function set.dim changes the dimension as wanted then
% changes the index to 0 and it calculates a new shiftValue.
        function set.dim(obj,newDim)%Setting a new dimension will change other properties(index and shiftValue)
           validateattributes(newDim,{'numeric'},{'scalar','positive','integer'});
            if newDim ~= obj.dim
               if newDim > length(obj.genVector)
                    disp('The dimension must be smaller than the length of the generating vector');
               else    
                   obj.privateDim=newDim;
                   obj.index=0;
                   obj.shiftValue = rand(1,newDim);
               end
           end
        end
        
%%
% set.genVector treats the input aswell. You can set genVector='default',
% and it will go back to the default generating vector. When you set a new
% vector, the function changes the index to 0 and set a new shiftValue.
% Also, if the current dimension is bigger than the new genVector, then the
% dimension will be changed to the length of the new genVector.
        function set.genVector(obj,newGenVector)
           if strcmp(newGenVector,'default')
               obj.privateGenVector = [1, 27808567, 14629093, 17951525, 873611, 20915065, 21016621, 23276977,...
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
           elseif length(newGenVector) ~= length(obj.genVector) || newGenVector ~= obj.genVector
              validateattributes(newGenVector,{'numeric'},{'vector','positive','integer'});
               obj.privateGenVector=newGenVector;
               obj.index=0;
               if length(newGenVector) < obj.dim
                    obj.dim= length(newGenVector);
               end
               obj.shiftValue = rand(1,obj.dim);
           end 

        end
        
        
    end
 
end
