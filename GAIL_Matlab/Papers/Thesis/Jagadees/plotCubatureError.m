
function plotCubatureError(dim, nvec, errCubMLE, ErrBound, fName, BernPolyOrder, ...
  ptransform, fullPath, visiblePlot, arbMean, s2, dsc)


if exist('visiblePlot','var') && visiblePlot==false
  hFigErr = figure('visible','off');
else
  hFigErr = figure();
end

% Bayesian Cubature Fourier Lattice
MATLABYellow = [0.9290, 0.6940, 0.1250];
set(hFigErr, 'units', 'inches', 'Position', [1 1 5.5 5.5])
%set(gca,'position',[0 0 1 1],'units','normalized')

%set(hFigErr, 'position',[0 0 1 1],'units','normalized')
%hold on

if exist('s2','var')
  if BernPolyOrder==4
    loglog(nvec,errCubMLE,'r.', ...
      nvec, dsc(1)*s2/s2(1), 'g-.', nvec, dsc, ':', ...
      nvec,abs(ErrBound),'b-', ...
      [nvec(1) nvec(end)],ErrBound(1)*[1 (nvec(1)/nvec(end))^2], '--') % color plots
    legend({'Actual Error', ...
      '\(\bar{s}\)', '\( (1-\frac{n}{\lambda_1})^{1/2} \)', ...
      'Error Bound', ...,
      '\(O(n^{-2})\)'}, ...
      'location','best')  % northeast
  else
    loglog(nvec,errCubMLE,'r.', ...
      nvec, dsc(1)*s2/s2(1), 'g-.', nvec, dsc, ':', ...
      nvec,abs(ErrBound),'b-', ...
      [nvec(1) nvec(end)],ErrBound(1)*[1 (nvec(1)/nvec(end))^1], '--') % color plots
    legend({'Actual Error', ...
      '\(\bar{s}\)', '\( (1-\frac{n}{\lambda_1})^{1/2} \)', ...
      'Error Bound', '\(O(n^{-1})\)'}, ...
      'location','best')  % northeast
  end
else
  if BernPolyOrder==4
    loglog(nvec,errCubMLE,'.', ...
      nvec,abs(ErrBound),'-.', ...
      [nvec(1) nvec(end)],errCubMLE(1)*[1 (nvec(1)/nvec(end))^4], '--', ...
      'color', MATLABYellow);
    legend({'Actual Error', ...
      'Error Bound', ...,
      '\(O(n^{-4})\)'}, ...
      'location','best')
  else
    loglog(nvec,errCubMLE,'.', ...
      nvec,abs(ErrBound),'-.', ...
      [nvec(1) nvec(end)],errCubMLE(1)*[1 (nvec(1)/nvec(end))^2], '--', ...
      'color', MATLABYellow);
    legend({'Actual Error', ...
      'Error Bound', '\(O(n^{-2})\)'}, ...
      'location','best')
  end
end

%title(sprintf('%s d=%d Bernoulli=%d, PeriodTx=%s', fName, dim, BernPolyOrder, ptransform))
if arbMean
  temp = '\(m \neq 0\)';
  title(sprintf('d=%d r=%d, Tx %s %s', dim, BernPolyOrder, ptransform, temp))
else
  title(sprintf('d=%d r=%d, Tx %s m=0', dim, BernPolyOrder, ptransform))
end

legend boxoff
if BernPolyOrder==4
  xlabel('Sample Size, \(n\)')
end
if dim==2
  ylabel('Error, \(|\mu - \hat{\mu}|\)')
end


% if (dim==3 && strcmp(fName, 'MVN')) || dim==4
%   xlabel('Sample Size, \(n\)')
% end
% if BernPolyOrder==2
%   ylabel('Error, \(|\mu - \hat{\mu}|\)')
% end


ylim([1E-18, 1E-1])
yticks([10^(-18) 10^(-15) 10^(-10) 10^(-5) 10^(-1)])  %(10.^[-18:5:-2])
xlim([2^10, 2^20])
%xticks([2^10 2^15 2^20])
xticks([10^4 10^5 10^6])
ytickangle(90)

%set(gca, 'xscale', 'log') % scale x-axis logarithmic
xticklabels({'\(10^{4}\)','\(10^{5}\)','\(10^{6}\)'})

plotFileName = sprintf('%s%s Error d_%d bernoulli_%d Period_%s.png', ...
  fullPath, fName, dim, BernPolyOrder, ptransform)  % let it print

set(gca,'LooseInset',get(gca,'TightInset'));
saveas(hFigErr, plotFileName)

end