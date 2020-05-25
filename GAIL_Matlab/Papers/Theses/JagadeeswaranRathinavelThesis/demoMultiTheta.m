
function demoMultiTheta()

nRep = 100;

muhatOneTheta = zeros(1,nRep);
const = [1E-4 1 1E4];
fun = @(x)sum(const.* sin(2*pi*x.^2), 2);
dim=3; absTol=1e-2;
exactInteg = fresnels(2)*sum(const)/2; % 'relTol',relTol, 
inputArgs = {'order',2, 'ptransform','none', ...
  'absTol',absTol, 'oneTheta',true,...
  'useGradient',false};
obj=cubBayesLattice_g(fun, dim, inputArgs{:});

[muhatOneTheta(nRep),outParamsOneTheta(nRep)]=compInteg(obj);
for i=1:nRep-1
  [muhatOneTheta(i),outParamsOneTheta(i)]=compInteg(obj);
end

mean(abs(muhatOneTheta-exactInteg))
mean([outParamsOneTheta.n])
mean([outParamsOneTheta.time])

muhatMultiTheta = zeros(1,nRep);
% outParamsMultiTheta(nRep) = {};
const = [1E-4 1 1E4];
fun = @(x)sum(const.* sin(2*pi*x.^2), 2);
dim=3; absTol=1e-2;
exactInteg = fresnels(2)*sum(const)/2; % 'relTol',relTol,
inputArgs = { 'order',2, 'ptransform','none', ...
  'absTol',absTol, 'oneTheta',false,...
  'useGradient',false};
obj=cubBayesLattice_g(fun, dim, inputArgs{:});
[muhatMultiTheta(nRep),outParamsMultiTheta(nRep)]=compInteg(obj);
for i=1:nRep-1
  [muhatMultiTheta(i),outParamsMultiTheta(i)]=compInteg(obj);
end

mean(abs(muhatMultiTheta-exactInteg))
mean([outParamsMultiTheta.n])
mean([outParamsMultiTheta.time])

outFileName = gail.save_mat('Paper_cubBayesLattice_g', 'MultiThetaDemo',...
    true, ...
    muhatMultiTheta, outParamsMultiTheta, ...
    muhatOneTheta, outParamsOneTheta);

MultiThetaDemoOut(outFileName)

fprintf('');
end

function MultiThetaDemoOut()

gail.InitializeDisplay

GAILPATH = GAILstart(0);
dirpath = strcat([GAILPATH,'OutputFiles',filesep], 'Paper_cubBayesLattice_g');

load(dataFileName);

avgAbsErrorOneTheta = mean(abs(muhatOneTheta-exactInteg));
avgSampleLenOneTheta = mean([outParamsOneTheta.n]);
avgTimeOneTheta = mean([outParamsOneTheta.time]);

avgAbsErrorMultiTheta =  mean(abs(muhatMultiTheta-exactInteg));
avgSampleLenMultiTheta = mean([outParamsMultiTheta.n]);
avgTimeMultiTheta = mean([outParamsMultiTheta.time]);


%% Output
fid = fopen([dirpath, filesep, 'MultiThetaOut.txt'], 'wt');

  fprintf(fid,'& \\multicolumn{2}{c}{\\text{Fresnel Sine Integral in} \\; d=3 } \\\\ \n');
  fprintf(fid,' \\hline \n');
  fprintf(fid,' \\text{Method} & \\text{OneTheta} & \\text{MultiTheta} \\\\ \n');
  fprintf(fid,' \\text{Absolute Error} & \\num{%1.5f} & \\num{%1.5f}   \\\\ \n', ...
     round([avgAbsErrorOneTheta; avgAbsErrorMultiTheta], 2, 'significant'));
  fprintf(fid,' n & \\num{%8.0f} & \\num{%8.0f} \\\\ \n', ...
     round([avgSampleLenOneTheta; avgSampleLenMultiTheta], 2, 'significant'));
  fprintf(fid,' \\text{Time (seconds)} & \\num{%7.4f} & \\num{%7.4f} \\\\ \n', ...
     round([avgTimeOneTheta; avgTimeMultiTheta], 2,'significant'));
  fprintf(fid,' \\hline');

fclose(fid);

end

