% generates the plots for the paper using given data or pre-computed data
function guaranteed_plots(varargin)

% close all
% if .mat file is not given then use the below files
if nargin < 1
  % myGAIL_path = GAILstart(0);
  % matFilePath = [myGAIL_path, filesep, ...
  %  'Papers\Thesis\Jagadees\Paper2018\figures\'];
  matFilePath = 'figures\';
  timestamp = '2018-Sep-6';
  
  filenames_ = {...
    'MVN_MLE_C2sin_d2_r2',...
    'MVN_full_C2sin_d2_r2',...
    'MVN_GCV_C2sin_d2_r2',...
    ...
    'Keister_MLE_C1sin_d4_r2',...
    'Keister_full_C1sin_d4_r2',...
    'Keister_GCV_C1sin_d4_r2',...
    ...
    'optPrice_MLE_Baker_d12_r1'...
    'optPrice_full_Baker_d12_r1'...
    'optPrice_GCV_Baker_d12_r1'...
    };
  % compose filenames by suffixing timestamp
  filenames = cellfun(@(a)['Guaranteed_plot_data_',a,'_',timestamp], ...
    filenames_, 'UniformOutput', false);
else
  matFilePath = varargin{1};
  filenames = varargin{2};
end

for filename=filenames
  generate_plot(matFilePath, filename{1});
end

end

function generate_plot(matFilePath, filename)

S = load([matFilePath filesep filename]);
figSavePath = matFilePath;
fig_size = [1 1 9 5];  % [1 1 9 6]
errTolVecText = arrayfun(@(x){sprintf('1e%d', x)}, S.log10ErrVec);

figHn = figure();
set(figHn, 'units', 'inches', 'Position', fig_size)
if strcmp(S.fName, 'MVN')
  errVecLimits = [1E-7, 1E1];
elseif strcmp(S.fName, 'Keister')
  errVecLimits = [1E-6, 1E1];
else
  errVecLimits = [1E-4, 1E1];
end
if isfield(S.outStructVec{1}(1), 'optParams')
  mvec = 2.^S.outStructVec{1}(1).optParams.mvec;
else
  mvec = 2.^S.outStructVec{1}(1).mvec;
end
nptsLimits = [(mvec(1)-1), (mvec(end)+1)];
nptsLimits(1) = 10^floor(log10(0.5*min(S.nptsVec(:))));
nptsLimits(2) = 10^ceil(log10(2*max(S.nptsVec(:))));

plot([1, 1], nptsLimits, 'r', 'LineWidth',1)
hold on
pointSize=30; %point size
pointShapes = {'<','s','>','d','^','v','o','p','h'};

for i=1:size(S.errVec,2)
  scatter(S.errVec(:,i),S.nptsVec(:,i),pointSize,log10(S.tolVec(:,i)),...
    pointShapes{i},'filled')
end

assert(max(S.nptsVec(:)) <= nptsLimits(2), ...
  sprintf('nume samples greater than max limit %d', nptsLimits(2)))

set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('\(\frac{\vert\mu-\widehat{\mu} \vert}{\varepsilon}\)','Interpreter','latex')
ylabel('Num. Samples')
c = colorbar('Direction','reverse', 'Ticks',S.log10ErrVec, ...
  'TickLabels',errTolVecText, 'TickLabelInterpreter','latex');
c.Label.Interpreter = 'latex';
c.Label.String = 'Error Tolerance, $\varepsilon$';
% axis tight; not required

axis([errVecLimits(1) errVecLimits(2) nptsLimits(1) nptsLimits(2)])
set(gca,'Xtick',(10.^(log10(errVecLimits(1)):4:log10(errVecLimits(2)))), ...
  'YTick',(10.^(floor(log10(nptsLimits(1))) :1:ceil(log10(nptsLimits(2))))))

if ~strcmp(S.testFunArg.sampling,'Lattice')
  S.testFunArg.varTx='';
end
figSavePathName = sprintf('%s_guaranteed_npts_%s_%s_d%d_r%d_%s.png', ...
  S.fName,S.stopCrit,S.testFunArg.varTx,...
  S.testFunArg.dim,S.testFunArg.order,S.timeStamp );
% gail.save_image(figHn, 'Paper_cubBayesLattice_g', figSavePathName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figH = figure();
set(figH, 'units', 'inches', 'Position', fig_size)

timeTicksLimits(1) = floor(log10(min(S.timeVec(:))));
timeTicksLimits(2) = ceil(log10(max(S.timeVec(:))));
plot([1, 1], 10.^timeTicksLimits, 'r', 'LineWidth',1)
hold on

% for i=1:size(S.errVec,2)
%   scatter(S.errVec(:,i),S.timeVec(:,i),pointSize,log10(S.tolVec(:,i)),...
%     pointShapes{i},'filled')
% end

set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('\({\vert\mu-\widehat{\mu} \vert}/{\varepsilon}\)', 'Interpreter', 'latex')
ylabel('Time (secs)')
c = colorbar('Direction','reverse', 'Ticks',S.log10ErrVec, ...
  'TickLabels',errTolVecText, 'TickLabelInterpreter','latex');
c.Label.Interpreter = 'latex';
c.Label.String = 'Error Tolerance, $\varepsilon$';


succeeded = find(S.exitflagVec(:,:) == 1);
  scatter(S.errVec(succeeded),S.timeVec(succeeded),pointSize,...
    log10(S.tolVec(succeeded)),'o','filled')
  missed = find(S.exitflagVec(:,:) ~= 1);
  scatter(S.errVec(missed),S.timeVec(missed),65,...
    log10(S.tolVec(missed)),'p')
% axis tight; not required
%axis([errVecLimits(1) errVecLimits(2) timeAxisLimits(1) timeAxisLimits(2) ])
%assert(max(timeVec(:)) <= timeLimits(2), sprintf('time val greater than max limit %d', timeLimits(2)))

axis([errVecLimits(1) errVecLimits(2) 10.^timeTicksLimits(1) 10.^timeTicksLimits(2) ])
timeTickInterval = ceil(diff(timeTicksLimits)/3);
set(gca,'Xtick',(10.^(log10(errVecLimits(1)):3:log10(errVecLimits(2)))), ...
  'YTick',(10.^(timeTicksLimits(1) :timeTickInterval: timeTicksLimits(2))))
%title(sprintf('%s d %d r %d %s %s', testFunArg.fName, ...
%              testFunArg.dim, testFunArg.order, testFunArg.varTx, mType));
figSavePathName = sprintf('%s_guaranteed_time_%s_%s_d%d_r%d_%s', ...
  S.fName,S.stopCrit,S.testFunArg.varTx,...
  S.testFunArg.dim,S.testFunArg.order,S.timeStamp );
gail.save_image(figH, 'Paper_cubBayesLattice_g', figSavePathName)

end
