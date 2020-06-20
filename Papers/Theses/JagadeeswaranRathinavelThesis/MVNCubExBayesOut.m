
function MVNCubExBayesOut(dataFileName)

% MVNCubExBayesOut: Prints LaTeX table of outputs from
% MVNCubatureExampleBayes

gail.InitializeDisplay

GAILPATH = GAILstart(0);
dirpath = strcat([GAILPATH,'OutputFiles',filesep], 'Paper_cubBayesLattice_g');

load(dataFileName);
%% Output
fid = fopen([dirpath, filesep, 'MVNBayesOut.txt'], 'wt');
for i = 1:length(avgAbsErrMC)
  if i == 1
  fprintf(fid,'& \\multicolumn{6}{c}{\\Sigma = \\bf{I}_d, \\ \\bf{b}=-\\bf{a}=(3.5,\\dots,3.5) } \\\\ \\hline \n');
  else
  fprintf(fid,'& \\multicolumn{6}{c}{\\Sigma = 0.4 \\bf{I}_d + 0.6\\bf{1}\\bf{1}^T , \\ \\bf{a}=(-\\infty,\\dots,-\\infty), \\ \\bf{b}=\\sqrt{d}(U_1,\\dots,U_d) } \\\\ \\hline \n');
  end
  fprintf(fid,' \\text{Method} & \\text{MC} & \\text{Lattice} & \\text{Sobol} & \\text{BayesLat} & \\text{BayesSobol}  \\\\\n');
  fprintf(fid,' \\text{Absolute Error} & \\num{%3.2e} & \\num{%3.2e} & \\num{%3.2e}  & \\num{%3.2e}  & \\num{%3.2e}  \\\\\n', ...
     round([avgAbsErrMC(i); avgAbsErrLat(i); avgAbsErrSob(i); avgAbsErrBayLat(i); avgAbsErrBaySob(i)], 2, 'significant'));
  fprintf(fid,' \\text{Tolerance Met} & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%%  \\\\\n', ...
     100*round([succMC(i); succLat(i); succSob(i); succBayLat(i); succBaySob(i)], 2, 'significant'));
  fprintf(fid,' n & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} \\\\\n', ...
     round([nSampleMC(i); nSampleLat(i); nSampleSob(i); nSampleBayLat(i); nSampleBaySob(i)], 2, 'significant'));
  fprintf(fid,' \\text{Time (seconds)} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f}  \n', ...
     round([timeMC(i); timeLat(i); timeSob(i); timeBayLat(i); timeBaySob(i)], 2,'significant'));
  if i == 1, fprintf(fid,' \\\\ \\\\'); end
end
fclose(fid);

end
