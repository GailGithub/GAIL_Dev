%This Matlab file may be used to download and install 'rbf_fd' into
% the current folder from where this code is run
%
% Step 1. Place this .m file where you want rbf_fd to go
%
% Step 2. Run this mfile

%% Download the package and change the directory
disp('The rbf_fd package is now downloading ...')

unzip('https://github.com/jagadeesr/rbf_fd/archive/v1_0.zip') %download ...and unzip
cd('rbf_fd-1_0') %get to the right subdirectory
if exist('rbfqr_v1_2', 'dir')~=7
    unzip('https://raw.github.com/jagadeesr/rbf_fd/master/rbfqr_mike_mccourt/rbfqr_v1_2.zip')
    cd('rbfqr_v1_2') % we need this package
    rbfsetup
    cd('..')
end

install_path

reply = input('\n Do you want to run script to test rbf_fd?\n y/n [n]:','s');
if any(strcmpi(reply,{'yes','y'}));
test_main
end
