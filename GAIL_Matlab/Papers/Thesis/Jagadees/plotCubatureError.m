
function plotCubatureError(dim, nvec, errCubMLE, ErrBound, fName, kernelOrder, ...
  ptransform, fullPath, visiblePlot, arbMean, scale, dsc)

k_order = kernelOrder/2;
k_order_str = sprintf('%1.1f',k_order);
k_order2_str = sprintf('%1.1f',k_order*2);
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

if exist('scale','var')
%   if kernelOrder==2
%     loglog(nvec,errCubMLE,'r.', ...
%       nvec, dsc(1)*scale/scale(1), 'g-.', nvec, dsc, ':', ...
%       nvec,abs(ErrBound),'b-', ...
%       [nvec(1) nvec(end)],ErrBound(1)*[1 (nvec(1)/nvec(end))^2], '--') % color plots
%     legend({'Actual Error', ...
%       '\(\bar{s}\)', '\( (1-\frac{n}{\lambda_1})^{1/2} \)', ...
%       'Error Bound', ...,
%       '\(O(n^{-2})\)'}, ...
%       'location','best')  % northeast
%   else
    loglog(nvec,errCubMLE,'r.', ...
      nvec, dsc(1)*scale/scale(1), 'g-.', nvec, dsc, ':', ...
      nvec,abs(ErrBound),'b-', ...
      [nvec(1) nvec(end)],ErrBound(1)*[1 (nvec(1)/nvec(end))^k_order], '--',...
      [nvec(1) nvec(end)],ErrBound(1)*[1 (nvec(1)/nvec(end))^(k_order*2)], '--' ...
    ) % color plots
    legend({'Actual Error', ...
      '\(\bar{s}\)', '\( (1-\frac{n}{\lambda_1})^{1/2} \)', ...
      'Error Bound', ['\(O(n^{-' k_order_str '})\)'], ...
      ['\(O(n^{-' k_order2_str '})\)']}, ...
      'location','best')  % northeast
%   end
else
%   if kernelOrder==2
    loglog(nvec,errCubMLE,'.', ...
      nvec,abs(ErrBound),'-.', ...
      [nvec(1) nvec(end)],errCubMLE(1)*[1 (nvec(1)/nvec(end))^k_order], '--', ...
      'color', MATLABYellow);
    legend({'Actual Error', ...
      'Error Bound', ['\(O(n^{-' k_order_str '})\)']}, ...
      'location','best')
%   else
%     loglog(nvec,errCubMLE,'.', ...
%       nvec,abs(ErrBound),'-.', ...
%       [nvec(1) nvec(end)],errCubMLE(1)*[1 (nvec(1)/nvec(end))^2], '--', ...
%       'color', MATLABYellow);
%     legend({'Actual Error', ...
%       'Error Bound', '\(O(n^{-2})\)'}, ...
%       'location','best')
%   end
end

%title(sprintf('%s d=%d Bernoulli=%d, PeriodTx=%s', fName, dim, kernelOrder, ptransform))
if arbMean
  temp = '\(m \neq 0\)';
  title(sprintf('d=%d order=%d, Tx %s %s', dim, kernelOrder, ptransform, temp))
else
  title(sprintf('d=%d r=%d, Tx %s m=0', dim, kernelOrder, ptransform))
end

legend boxoff
if kernelOrder==4
  xlabel('Sample Size, \(n\)')
end
if dim==2
  ylabel('Error, \(|\mu - \hat{\mu}|\)')
end


% if (dim==3 && strcmp(fName, 'MVN')) || dim==4
%   xlabel('Sample Size, \(n\)')
% end
% if kernelOrder==2
%   ylabel('Error, \(|\mu - \hat{\mu}|\)')
% end


ylim([1E-18, 1E-0])
yticks(10.^[-15 -10 -5 -1])  %(10.^[-18:5:-2])
ytickangle(90)

xlim([10^3, 2^20])
%xticks([2^10 2^15 2^20])
xticks(10.^[3 4 5 6])

%set(gca, 'xscale', 'log') % scale x-axis logarithmic
xticklabels({'\(10^{3}\)','\(10^{4}\)','\(10^{5}\)','\(10^{6}\)'})
grid on

plotFileName = sprintf('%s%s Error d_%d order_%d Period_%s.png', ...
  fullPath, fName, dim, kernelOrder, ptransform)  % let it print

set(gca,'LooseInset',get(gca,'TightInset'));
saveas(hFigErr, plotFileName)

end