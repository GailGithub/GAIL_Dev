% KeisterCubExBayesOut: Prints LaTeX table of outputs from
% KeisterCubatureExampleBayes
function KeisterCubExBayesOut(dataFileName)

gail.InitializeDisplay

GAILPATH = GAILstart(0);
dirpath = strcat([GAILPATH,'OutputFiles',filesep], 'Paper_cubBayesLattice_g');

load(dataFileName);

%% Output
fid = fopen([dirpath, filesep, 'KeisterBayesOut.txt'], 'wt');
for i = 1:length(avgAbsErrMC)
  fprintf(fid,'& \\multicolumn{4}{c}{d = %3d,\\ \\varepsilon = %6.3f} \\\\ \n', dvec(i), abstol(i));
  fprintf(fid,' \\hline \n');
  fprintf(fid,' \\text{Method} & \\text{MC} & \\text{Lattice} & \\text{Sobol} & \\text{BayesLat} & \\text{BayesSobol}  \\\\ \n');
  fprintf(fid,' \\text{Absolute Error} & \\num{%1.6f} & \\num{%1.6f} & \\num{%1.6f}  & \\num{%1.6f}  & \\num{%1.6f}  \\\\ \n', ...
     round([avgAbsErrMC(i); avgAbsErrLat(i); avgAbsErrSob(i); avgAbsErrBayLat(i); avgAbsErrBaySob(i)], 2, 'significant'));
  fprintf(fid,' \\text{Tolerance Met} & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%%  \\\\\n', ...
     100*round([succMC(i); succLat(i); succSob(i); succBayLat(i); succBaySob(i)], 2, 'significant'));
  fprintf(fid,' n & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} \\\\\n', ...
     round([nSampleMC(i); nSampleLat(i); nSampleSob(i); nSampleBayLat(i); nSampleBaySob(i)], 2, 'significant'));
  fprintf(fid,' \\text{Time (seconds)} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} \\\\ \n', ...
     round([timeMC(i); timeLat(i); timeSob(i); timeBayLat(i); timeBaySob(i)], 2,'significant'));
  fprintf(fid,' \\\\ \n');
end
fprintf(fid,' \\hline \n');
fclose(fid);

