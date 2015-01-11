function GAIL_Install()
% GAIL_INSTALL   Install GAIL. Add GAIL paths to MATLAB search path.
[GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION] = GAILstart;
fprintf('\nWelcome to Guaranteed Automatic Integration Library (GAIL).\nYou are installing GAIL v%s.\n\n', GAILVERSION);
if MATLABVERSION < 7,
  error('This version is only supported on Matlab 7.x and above.');
end
gailp=genpath(GAILPATH);% Generate strings of paths to GAIL subdirectories
warninfo = warning('query','MATLAB:rmpath:DirNotFound');
warning('off',warninfo.identifier);
rmpath(gailp);% Remove path from MATLAB search path
warning(warninfo.state, warninfo.identifier);
addpath(gailp);           % Add GAIL directories and subdirectories
savepath;                 % Save the changes
warninfo = warning('query','MATLAB:doc:DocNotInstalled');
warning('off', warninfo.identifier);
GAIL_Publish; 
builddocsearchdb(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html'));
warning(warninfo.state, warninfo.identifier);
fprintf('\nYou can go to help documentation ---> Supplemental Software to learn how to use GAIL.\n');
%addpath(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html',PATHNAMESEPARATOR,'helpsearch'));
fprintf('\nGAIL version %s has been installed successfully.\n\n', GAILVERSION);
end

function GAIL_Publish()
% GAIL_PUBLISH  To generate html files in the GAIL subdirectory Documentation
oldStatus = get(0,'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off')
publish('GAIL');
publish('funclist');
publish('help_funappx_g');
publish('help_integral_g');
publish('help_meanMC_g');
publish('help_meanMCBer_g');
publish('help_cubMC_g');
publish('help_funmin_g');
publish('help_cubLattice_g');
publish('help_cubSobol_g');
publish('help_install');
set(0, 'DefaultFigureVisible', oldStatus)
end
