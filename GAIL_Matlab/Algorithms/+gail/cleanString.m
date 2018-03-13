function outString = cleanString(inString,EPref)
% This function takes strings of MATLAB formatted numbers and replaces
% certain characters by others that look better for a table
outString = strrep(inString,'E','\text{E}');
outString = strrep(outString,'-0','{-}');
outString = strrep(outString,'+0','');
if nargin < 2
   EPref = 'QQ';
end
outString = strrep(outString,[EPref '\text{E}'],'E'); %allow an E in the text