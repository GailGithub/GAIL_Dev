function guaranteed_plots(varargin)

if nargin < 1
matFilePath = 'D:\Dropbox\fjhickernellGithub\GAIL_Dev-BayesianCubature\GAIL_Matlab\Papers\Thesis\Jagadees\Paper2018\figures\';

% filenames = {...
%   'Guaranteed_plot_data_MVN_MLE_C2sin_d3_r4_2018-Aug-19', ...
%   'Guaranteed_plot_data_MVN_full_C2sin_d3_r4_2018-Aug-19', ...
%   'Guaranteed_plot_data_MVN_GCV_C2sin_d3_r4_2018-Aug-19', ...
%   ...
%   'Guaranteed_plot_data_Keister_MLE_C1_d4_r4_2018-Aug-19',...
%   'Guaranteed_plot_data_Keister_full_C1_d4_r4_2018-Aug-19',...
%   'Guaranteed_plot_data_Keister_GCV_C1_d4_r4_2018-Aug-19',...
%   ...
%   'Guaranteed_plot_data_optPrice_MLE_Baker_d12_r2_2018-Aug-19'...
%   'Guaranteed_plot_data_optPrice_full_Baker_d12_r2_2018-Aug-19'...
%   'Guaranteed_plot_data_optPrice_GCV_Baker_d12_r2_2018-Aug-19'...
%   };

filenames = {...
  'Guaranteed_plot_data_MVN_MLE_C2sin_d2_r4_2018-Aug-26',...
  'Guaranteed_plot_data_MVN_full_C2sin_d2_r4_2018-Aug-26',...
  'Guaranteed_plot_data_MVN_GCV_C2sin_d2_r4_2018-Aug-26',...
  ...
  'Guaranteed_plot_data_Keister_MLE_C1sin_d4_r2_2018-Aug-26',...
  'Guaranteed_plot_data_Keister_full_C1sin_d4_r2_2018-Aug-26',...
  'Guaranteed_plot_data_Keister_GCV_C1sin_d4_r2_2018-Aug-26',...
  };
else
  matFilePath = varargin{1};
  filenames = varargin{2};
end

for filename=filenames
  filename=filename{1};
  generate_plot(matFilePath, filename);
end

end

function generate_plot(matFilePath, filename)

S = load([matFilePath filename]);
figSavePath = matFilePath;

errTolVecText = arrayfun(@(x){sprintf('1e%d', x)}, S.log10ErrVec);

figHn = figure();
set(figHn, 'units', 'inches', 'Position', [1 1 9 6])
errVecLimits = [1E-7, 1E1];
mvec = 2.^S.outStructVec{1}(1).mvec;
nptsLimits = [(mvec(1)-1), (mvec(end)+1)];
nptsLimits(1) = 10^floor(log10(0.5*min(S.nptsVec(:))));
nptsLimits(2) = 10^ceil(log10(2*max(S.nptsVec(:))));

plot([1, 1], nptsLimits, 'r', 'LineWidth',1)
hold on
pointSize=30; %point size
pointShapes = {'o','s','d','^','v','<','>','p','h'};

for i=1:size(S.errVec,2)
  scatter(S.errVec(:,i),S.nptsVec(:,i),pointSize,log10(S.tolVec(:,i)),pointShapes{i},'filled')
end

assert(max(S.nptsVec(:)) <= nptsLimits(2), ...
  sprintf('nume samples greater than max limit %d', nptsLimits(2)))

set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('\(\frac{\vert\mu-\widehat{\mu} \vert}{\varepsilon}\)')
ylabel('Num. Samples')
c = colorbar('Direction','reverse', 'Ticks',S.log10ErrVec, ...
  'TickLabels',errTolVecText, 'TickLabelInterpreter','tex');
c.Label.Interpreter = 'latex';
c.Label.String = 'Error Tolerance, $\varepsilon$';
% axis tight; not required

axis([errVecLimits(1) errVecLimits(2) nptsLimits(1) nptsLimits(2)])
set(gca,'Xtick',(10.^(log10(errVecLimits(1)):4:log10(errVecLimits(2)))), ...
  'YTick',(10.^(floor(log10(nptsLimits(1))) :1:ceil(log10(nptsLimits(2))))))
% if testFunArg.arbMean
%   mType = '\(m \neq 0\)'; % arb mean
% else
%   mType = '\(m = 0\)'; % zero mean
% end
%   title(sprintf('%s d %d r %d %s %s', testFunArg.fName, ...
%                 testFunArg.dim, testFunArg.order, testFunArg.varTx, mType));

figSavePathName = sprintf('%s%s_guaranteed_npts_%s_%s_d%d_r%d_%s.png', ...
  figSavePath, S.fName,S.stopCrit,S.testFunArg.varTx,...
  S.testFunArg.dim,S.testFunArg.order,S.timeStamp );
% saveas(figHn, figSavePathName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figH = figure();
set(figH, 'units', 'inches', 'Position', [1 1 9 6])
%timeLimits = [1E-3, 1E-2];
timeAxisLimits(1) = min(S.timeVec(:))/2;
timeAxisLimits(2) = max(S.timeVec(:))*2;
plot([1, 1], timeAxisLimits, 'r', 'LineWidth',1)
hold on

for i=1:size(S.errVec,2)
  scatter(S.errVec(:,i),S.timeVec(:,i),pointSize,log10(S.tolVec(:,i)),...
    pointShapes{i},'filled')
end

set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('\(\frac{\vert\mu-\widehat{\mu} \vert}{\varepsilon}\)')
ylabel('Time (secs)')
c = colorbar('Direction','reverse', 'Ticks',S.log10ErrVec, ...
  'TickLabels',errTolVecText, 'TickLabelInterpreter','tex');
c.Label.Interpreter = 'latex';
c.Label.String = 'Error Tolerance, $\varepsilon$';
% axis tight; not required

%assert(max(timeVec(:)) <= timeLimits(2), sprintf('time val greater than max limit %d', timeLimits(2)))
axis([errVecLimits(1) errVecLimits(2) timeAxisLimits(1) timeAxisLimits(2) ])

timeTicksLimits(1) = ceil(log10(min(S.timeVec(:)))); %floor
timeTicksLimits(2) = ceil(log10(max(S.timeVec(:))));
set(gca,'Xtick',(10.^(log10(errVecLimits(1)):3:log10(errVecLimits(2)))), ...
  'YTick',(10.^(timeTicksLimits(1) :1: timeTicksLimits(2))))
%title(sprintf('%s d %d r %d %s %s', testFunArg.fName, ...
%              testFunArg.dim, testFunArg.order, testFunArg.varTx, mType));  
figSavePathName = sprintf('%s%s_guaranteed_time_%s_%s_d%d_r%d_%s.png', ...
  figSavePath, S.fName,S.stopCrit,S.testFunArg.varTx,...
  S.testFunArg.dim,S.testFunArg.order,S.timeStamp );
saveas(figH, figSavePathName)

fprintf('done\n')
end