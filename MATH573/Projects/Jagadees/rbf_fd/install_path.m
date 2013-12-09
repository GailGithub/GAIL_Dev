    
rbfdversion = 1.0;
fprintf('\nWelcome to rbf_fd version %2.1f \n', rbfdversion);
Friend = computer;

if isunix % for Mac or unix
    PATHNAMESEPARATOR = '/';
    RBFFDPATH=[pwd, PATHNAMESEPARATOR];
elseif strcmp(Friend(1:2),'PC') % for pc
    PATHNAMESEPARATOR = '\';
    RBFFDPATH = [pwd, PATHNAMESEPARATOR];
else
    error('I don''t recognize this computer.')
end


V = version;
wh = strfind(V,'.');
wh = wh(2)-1;
MATLABVERSION = str2double(V(1:wh));
if MATLABVERSION < 7,
  error('This version is only supported on Matlab 7.x and above.');
else
  rbffdp=genpath(RBFFDPATH);% adding all subdirectories 
end
addpath(rbffdp);% Add directories and subdirectories
savepath; % Save the changes
fprintf('\n rbf_fd version %g has been installed successfully.\n\n', rbfdversion);

