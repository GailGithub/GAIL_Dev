%GAILBUILD To build doc search database in MATLAB help
[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
helppath = strcat(GAILPATH,'Documentation' , PATHNAMESEPARATOR, 'html');
builddocsearchdb(helppath)

