function saveMCQMC2012peakyfundir(fname,xsample)
GAILPATH = GAILstart(0);
outputfolder =  [GAILPATH,'OutputFiles',filesep,'MCQMC2012PaperOutput'];
if exist(outputfolder) ~= 7,
  mkdir(outputfolder);
end
fullfilename = strcat(outputfolder, filesep', fname);
save(fullfilename,'xsample')
end