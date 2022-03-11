% KeisterCubExWileyOut: Prints LaTeX table of outputs from KeisterCubatureExampleWiley
function KeisterCubExJorsOut(dataFileName)

gail.InitializeDisplay

GAILPATH = GAILstart(0);
dirpath = strcat([GAILPATH,'OutputFiles',filesep], 'GAIL_JORS2021_Output');

load(dataFileName);

%% Output
fid = fopen([dirpath, filesep, 'KeisterJorsOut.txt'], 'wt');
for i = 1:length(avgAbsErrMC)
  fprintf(fid,'& \\multicolumn{5}{c}{d = %3d,\\ \\varepsilon = %6.3f} \\\\ \\hline \n', dvec(i), abstol(i));
  fprintf(fid,' \\text{Method} & \\text{MC} & \\text{Lattice} & \\text{Sobol} & \\text{Bayes Lattice} & \\text{Bayes Net}  \\\\\n');
  fprintf(fid,' \\text{Absolute Error} & \\num{%7.1e} & \\num{%7.1e} & \\num{%7.1e}  & \\num{%7.1e}  & \\num{%7.1e}  \\\\\n', ...
     round([avgAbsErrMC(i); avgAbsErrLat(i); avgAbsErrSob(i); avgAbsErrBay(i); avgAbsErrBayNet(i)], 2, 'significant'));
  fprintf(fid,' \\text{Tolerance Met} & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%%  & \\num{%5.0f} \\%%  \\\\\n', ...
     100*round([succMC(i); succLat(i); succSob(i); succBay(i); succBayNet(i)], 2, 'significant'));
  fprintf(fid,' n & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f}  & \\num{%8.0f}  \\\\\n', ...
     round([nSampleMC(i); nSampleLat(i); nSampleSob(i); nSampleBay(i); nSampleBayNet(i)], 2, 'significant'));
  fprintf(fid,' \\text{Time (seconds)} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} \n', ...
     round([timeMC(i); timeLat(i); timeSob(i); timeBay(i); timeBayNet(i)], 2,'significant'));
  if i == 1, fprintf(fid,' \\\\ '); end
end
fclose(fid);

