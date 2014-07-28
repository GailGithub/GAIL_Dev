%GAILPUBLISH  To generate html files in the GAIL subdirectory Documentation
publish('GAIL');
publish('funclist');
publish('help_funappx_g');
publish('help_integral_g');
publish('help_meanMC_g');
publish('help_cubMC_g');
[GAILPATH,~,PATHNAMESEPARATOR,~] = GAILstart(0);
builddocsearchdb(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html'));