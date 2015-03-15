function saveMCQMC2012peakyfundir(fname,xsample)
[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
outputfolder =  [GAILPATH,'OutputFiles',PATHNAMESEPARATOR,'MCQMC2012PaperOutput'];
if exist(outputfolder) ~= 7,
  mkdir(outputfolder);
end
fullfilename = strcat(outputfolder, PATHNAMESEPARATOR', fname);
save(fullfilename,'xsample')
end