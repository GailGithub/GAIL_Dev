%This mfile may be used to download and install XuanZhouQeQ into 
%   the location you choose
%
%   Step 1.  Place this .m file in where you want XuanZhouQeQ to go
%
%   Step 2.  Run this mfile

%% Download the package and change the directory
disp('Downloading and installing XuanZhouQeq...');
unzip(['https://docs.google.com/uc?id=0Bz5Mny5Fq8w0VUFOandwNjlILUE&export=download']) %download and unzip
%https://drive.google.com/file/d/0Bz5Mny5Fq8w0VUFOandwNjlILUE/edit?usp=sharing
cd('XuanZhouQeq') %get to the right subdirectory

%% Install XuanZhouQeq
clear all; close all; 
Friend = computer;
if isunix % for Mac or unix
    PATHNAMESEPARATOR = '/';
    XuanZhouQeqPATH=[fileparts(which('XuanZhouDownloadInstall')), PATHNAMESEPARATOR];
elseif strcmp(Friend(1:2),'PC') % for pc
    PATHNAMESEPARATOR = '\';
    XuanZhouQeqPATH = [fileparts(which('XuanZhouDownloadInstall')), PATHNAMESEPARATOR];
else
    error('I don''t recognize this computer.')
end

XuanZhouQeqp=genpath(XuanZhouQeqPATH);% adding all subdirectories 
addpath(XuanZhouQeqp);% Add XuanZhouQeq directories and subdirectories
savepath; % Save the changes

%% Run a quick test
XuanZhouQeqWorkout;
%   Then you should be ready to use XuanZhouQeq