% GAILSTART  Initialize all the GAIL paths and parameters
function [GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION] = GAILstart(isverbose)

if nargin < 1
    isverbose =  true;
end % print the variable names and values

GAILVERSION = 1.3;

Friend = computer;

if isunix % for Mac or unix
    PATHNAMESEPARATOR = '/';
    GAILPATH=[fileparts(which('GAILstart')), PATHNAMESEPARATOR];
elseif strcmp(Friend(1:2),'PC') % for pc
    PATHNAMESEPARATOR = '\';
    GAILPATH = [fileparts(which('GAILstart')), PATHNAMESEPARATOR];
else
    error('I don''t recognize this computer.')
end

V = version;
wh = strfind(V,'.');
wh = wh(2)-1;
MATLABVERSION = str2double(V(1:wh));
if(isverbose)
    fprintf('Setting Output:\n');
    fprintf('   MATLABVERSION = %g\n',        MATLABVERSION)
    fprintf('   GAILVERSION = %g\n',          GAILVERSION)
    fprintf('   GAILPATH = %s\n',             GAILPATH)
    fprintf('   PATHNAMESEPARATOR = "%s"\n',  PATHNAMESEPARATOR)
end