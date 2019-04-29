% KeisterCubExWileyOut: Prints LaTeX table of outputs from KeisterCubatureExampleWiley
function KeisterCubExPearcOut(dataFileName)

gail.InitializeDisplay

GAILPATH = GAILstart(0);
dirpath = strcat([GAILPATH,'OutputFiles',filesep], 'GAIL_PEARC2019_Output');

load(dataFileName);

%% Output
fid = fopen([dirpath, filesep, 'KeisterPearcOut.txt'], 'wt');
for i = 1:length(avgAbsErrMC)
  fprintf(fid,'& \\multicolumn{4}{c}{d = %3d,\\ \\varepsilon = %6.3f} \\\\ \\hline \n', dvec(i), abstol(i));
  fprintf(fid,' \\text{Method} & \\text{MC} & \\text{Lattice} & \\text{Sobol} & \\text{Bayes}  \\\\\n');
  fprintf(fid,' \\text{Absolute Error} & \\num{%7.5f} & \\num{%7.5f} & \\num{%7.5f}  & \\num{%7.5f}  \\\\\n', ...
     round([avgAbsErrMC(i); avgAbsErrLat(i); avgAbsErrSob(i); avgAbsErrBay(i)], 2, 'significant'));
  fprintf(fid,' \\text{Tolerance Met} & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%%  \\\\\n', ...
     100*round([succMC(i); succLat(i); succSob(i); succBay(i)], 2, 'significant'));
  fprintf(fid,' n & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f}  \\\\\n', ...
     round([nSampleMC(i); nSampleLat(i); nSampleSob(i); nSampleBay(i)], 2, 'significant'));
  fprintf(fid,' \\text{Time (seconds)} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f}  \n', ...
     round([timeMC(i); timeLat(i); timeSob(i); timeBay(i)], 2,'significant'));
  if i == 1, fprintf(fid,' \\\\ \\\\'); end
end
fclose(fid);

