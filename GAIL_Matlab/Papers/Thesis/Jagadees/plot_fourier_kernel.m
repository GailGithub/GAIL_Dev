% plots the shift invariant fourier kernel for demo purposes
function plot_fourier_kernel()

n = 100;

% using uniform points
l = linspace(0,1,n);
[X, Y] = meshgrid(l);

shape_param = [0.1 0.9];
order = [2 4];

for sh = shape_param
  for r = order
    hFig = figure('visible','off');
    set(hFig, 'units', 'inches', 'Position', [4 4 6.5 5.5])
    Z = kernel([X(:) Y(:)], r, sh);
    meshc(X, Y, reshape(Z, [n, n]));
    
    title(sprintf('r=%d shape=%0.2f', r, sh))
    figSavePathName = sprintf('fourier_kernel r_%d shape_%dby100.png', r, 100*sh)
    saveas(hFig, figSavePathName)
  end
end

end


% bernoulli polynomial based kernel
function [K] = kernel(xun,order,a)

constMult = -(-1)^(order/2)*((2*pi)^order)/factorial(order);
if order == 2
  bernPloy = @(x)(-x.*(1-x) + 1/6);
elseif order == 4
  bernPloy = @(x)( ( (x.*(1-x)).^2 ) - 1/30);
else
  error('Bernoulli order not implemented !');
end

K = prod(1 + (a)*constMult*bernPloy(xun),2);

end
