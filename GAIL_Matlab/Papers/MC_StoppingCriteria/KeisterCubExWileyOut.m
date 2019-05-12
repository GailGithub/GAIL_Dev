% KeisterCubExWileyOut: Prints LaTeX table of outputs from KeisterCubatureExampleWiley
function KeisterCubExWileyOut(dataFileName)

gail.InitializeDisplay

GAILPATH = GAILstart(0);
dirpath = strcat([GAILPATH,'OutputFiles',filesep], 'MC_StoppingCriteriaOutput');

load(dataFileName);

%% Output
fid = fopen([dirpath, filesep, 'KeisterOut.txt'], 'wt');
fprintf(fid,'& \\multicolumn{3}{c}{d = 3,\\ \\varepsilon = %6.3f} & \\multicolumn{3}{c}{d = 8,\\ \\varepsilon = %6.2f} \\\\\n',abstol);
fprintf(fid,' \\text{Method} & \\texttt{cubMC\\_g}  & \\texttt{cubLattice\\_g} & \\texttt{cubSobol\\_g} & \\texttt{cubMC\\_g} & \\texttt{cubLattice\\_g} & \\texttt{cubSobol\\_g}  \\\\\n');
fprintf(fid,' \\text{Absolute Error} & \\num{%7.5f} & \\num{%7.5f} & \\num{%7.5f} & \\num{%7.5f} & \\num{%7.5f} & \\num{%7.5f}  \\\\\n', ...
   round([avgAbsErrMC; avgAbsErrLat; avgAbsErrSob],2,'significant'));
fprintf(fid,' \\text{Tolerance Met} & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%% & \\num{%5.0f} \\%%  \\\\\n', ...
   100*round([succMC; succLat; succSob],2,'significant'));
fprintf(fid,' n & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} & \\num{%8.0f} \\\\\n', ...
   round([nSampleMC; nSampleLat; nSampleSob],2,'significant'));
fprintf(fid,' \\text{Time (seconds)} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f} & \\num{%7.4f}  \n', ...
   round([timeMC; timeLat; timeSob],2,'significant'));
fclose(fid);

