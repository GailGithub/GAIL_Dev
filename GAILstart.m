% GAILSTART  Initialize all the GAIL paths and system parameters.
function [GAILPATH,GAILVERSION,MATLABVERSION] = GAILstart(isverbose)

if nargin < 1
    isverbose =  true;
end % print the variable names and values

GAILVERSION = '2.3.2';
GAILPATH=[fileparts(which('GAILstart')),filesep];
V = version;
wh = strfind(V,'.');
wh = wh(2)-1;
MATLABVERSION = str2double(V(1:wh));
if(isverbose)
    fprintf('Setting Output:\n');
    fprintf('   MATLABVERSION = %g\n',        MATLABVERSION)
    fprintf('   GAILVERSION = %s\n',          GAILVERSION)
    fprintf('   GAILPATH = %s\n',             GAILPATH)
end
