
function plotCubatureError(dim, nvec, errCubMLE, ErrBound, fName, BernPolyOrder, ptransform, fullPath, visiblePlot)
if exist('visiblePlot','var') && visiblePlot==false
    hFigErr = figure('visible','off');
else
    hFigErr = figure(); %21
end
    MATLABYellow = [0.9290, 0.6940, 0.1250];
    set(hFigErr, 'units', 'inches', 'Position', [4 4 5.5 5.5])
    %set(gca,'position',[0 0 1 1],'units','normalized')
    
    %set(hFigErr, 'position',[0 0 1 1],'units','normalized')
    %hold on
    
    if true
        if BernPolyOrder==4
            loglog(nvec,errCubMLE,'.', ...
                nvec,abs(ErrBound),'-.', ...
                [nvec(1) nvec(end)],errCubMLE(1)*[1 (nvec(1)/nvec(end))^4], '--', ...
                'color', MATLABYellow);
            legend({'Actual Error', ...
                'Error Bound', ...,
                '\(O(n^{-4})\)'}, ...
                'location','southwest')
        else
            loglog(nvec,errCubMLE,'.', ...
                nvec,abs(ErrBound),'-.', ...
                [nvec(1) nvec(end)],errCubMLE(1)*[1 (nvec(1)/nvec(end))^2], '--', ...
                'color', MATLABYellow);
            legend({'Actual Error', ...
                'Error Bound', '\(O(n^{-2})\)'}, ...
                'location','southwest')
            % Bayesian Cubature Fourier Lattice 
        end
    else
        
        if BernPolyOrder==4
            loglog(nvec,errCubMLE,'.', ...
                nvec,abs(ErrBound),'-.', ...
                [nvec(1) nvec(end)],errCubMLE(1)*[1 (nvec(1)/nvec(end))^4], '--', ...
                [nvec(1) nvec(end)],errCubMLE(1)*[1 (nvec(1)/nvec(end))^2], ':', ...
                'color', MATLABYellow);
            legend({'Actual Error', ...
                'Error Bound', ...,
                '\(O(n^{-4})\)', '\(O(n^{-2})\)'}, ...
                'location','southwest')
        else
            loglog(nvec,errCubMLE,'.', ...
                nvec,abs(ErrBound),'-.', ...
                [nvec(1) nvec(end)],errCubMLE(1)*[1 (nvec(1)/nvec(end))^2], '--', ...
                [nvec(1) nvec(end)],errCubMLE(1)*[1 (nvec(1)/nvec(end))^1], ':', ...
                'color', MATLABYellow);
            legend({'Actual Error', ...
                'Error Bound', '\(O(n^{-2})\)', '\(O(n^{-1})\)'}, ...
                'location','southwest')
            % Bayesian Cubature Fourier Lattice 
        end
    end
    
    %title(sprintf('%s d=%d Bernoulli=%d, PeriodTx=%s', fName, dim, BernPolyOrder, ptransform))
    title(sprintf('d=%d r=%d, Tx=%s', dim, BernPolyOrder, ptransform))
    
    legend boxoff
    if BernPolyOrder==4
        xlabel('Sample Size, \(n\)')
    end
    if dim==2
        ylabel('Error, \(|\mu - \hat{\mu}|\)')
    end
    
    ylim([1E-18, 1E-1])
    yticks([10^(-18) 10^(-15) 10^(-10) 10^(-5) 10^(-1)])  %(10.^[-18:5:-2])
    xlim([2^10, 2^23])
    xticks([2^10 2^15 2^20 2^23])
    ytickangle(90)
    
    
    %set(gca, 'xscale', 'log') % scale x-axis logarithmic
    xticklabels({'\(2^{10}\)','\(2^{15}\)','\(2^{20}\)','\(2^{23}\)'})
    
    set(gca,'LooseInset',get(gca,'TightInset'));
    saveas(hFigErr, sprintf('%s%s Error d_%d bernoulli_%d Period_%s.png', ...
        fullPath, fName, dim, BernPolyOrder, ptransform))

end