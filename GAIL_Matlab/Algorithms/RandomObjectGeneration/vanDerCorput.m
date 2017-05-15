function out = vanDerCorput(varargin)
%% vanDerCorput
% This function permit you to enter 2 indexes, what will return all Van Der Corput
% numbers between these indexes. Or you can enter one index
% and it will be the number of Van Der Corput elements returned starting from 0.
%   Example:
%       >>vanDerCorput(1,5)
%       ans =
%           0.5000
%           0.2500
%           0.7500
%           0.1250
%           0.6250
%
%Author: Maurício Darós Andrade
%The idea was mainly taken from http://rosettacode.org/wiki/Van_der_Corput_sequence

    if nargin >= 2
        start = varargin{1};%The index of the first element (0th element = 0; 1st element = 0.5 ...)
        final = varargin{2};%The index of the last element
    else
        start = 0;
        final = varargin{1}-1;%here the input is the amount of elements wanted
    end
    validateattributes(start,{'numeric'},{'scalar','nonnegative','integer'})
    validateattributes(final,{'numeric'},{'scalar','nonnegative','integer'})
    assert(start<=final,'The last index must be greater or equal to the first');
    
    binary = dec2bin(start:final)-'0';%binary is a matrix in which each line ...
    %is a number from start to final. Each element of a line is a digit of
    %the number in binary. I.e, for start=1,final=3 then binary= 0 0 1
    %                                                            0 1 0
    %                                                            0 1 1
    %                                                            1 0 0
    %The -'0' is nescessary to transform the vector in number since the
    %dec2bin transform decimal numbers into binary strings. '0' in ASCII is
    %48 while '1' is 49. Therefore '0' get to be 0 and '1'-'0' get to be 1.
    width = size(binary,2); %width is the number of columns 
    v = (1:width) - width - 1; %v is a vector from -width untill -1
    out = binary*(2.^v'); %2.^v' are the weights for the binary digits.
end


