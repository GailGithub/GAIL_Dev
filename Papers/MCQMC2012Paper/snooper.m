function y=snooper(x,info)
%This records all the values of x it is asked for in a file named what you
%like (.txt, .mat, or other suffix must be specially typed)
GAILPATH = GAILstart(0);
outputfolder =  [GAILPATH,'OutputFiles',filesep,'MCQMC2012PaperOutput'];
fullfilename = strcat(outputfolder, filesep', info.filename);
load(fullfilename,'xsample')
xsample=[xsample; x(:)];
save(fullfilename,'xsample')
if(info.RegFunc(x)==0)
    y=zeros(size(x));
else
y=info.RegFunc(x); %This is the function that you are trying to fit peaks into
end


