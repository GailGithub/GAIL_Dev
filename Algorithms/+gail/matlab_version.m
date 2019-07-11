function MATLABVERSION = matlab_version
V = version;
wh = strfind(V,'.');
wh = wh(2)-1;
MATLABVERSION = str2double(V(1:wh));