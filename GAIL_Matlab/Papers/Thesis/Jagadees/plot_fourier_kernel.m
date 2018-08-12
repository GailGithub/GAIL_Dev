% plots the shift invariant fourier kernel for demo purposes
function plot_fourier_kernel()

n = 512;
shape_param = [0.2 0.8];
order = [2 4];
dim = 1;

% using uniform points
l = linspace(0,1,n);

if dim==1
  xpts = l';
  Z = zeros(size(xpts, 1), 4);
  hFig = figure('visible','on');
  set(hFig, 'units', 'inches', 'Position', [4 4 6.5 5.5])
  hold on
  i=1;
  leg_text = cell(4,1);
  for r = order
    for sh = shape_param
      Z(:,i) = kernel(xpts, r, sh);
      leg_text{i} = sprintf('$r=%d,\\theta=%1.1f$', r, sh);
      i=i+1;
    end
  end

  plot(xpts,Z)
  xlabel('\(x\)')
  ylabel('\(C_{\theta}(x,0.3)\)')
  legend(leg_text, 'Interpreter','latex','location','best')
  axis tight

  figSavePathName = sprintf('fourier_kernel dim_1.png');
  saveas(hFig, figSavePathName)
end

if dim==2
  [X, Y] = meshgrid(l);
  
  for sh = shape_param
    for r = order
      hFig = figure('visible','off');
      set(hFig, 'units', 'inches', 'Position', [4 4 6.5 5.5])
      Z = kernel([X(:) Y(:)], r, sh);
      meshc(X, Y, reshape(Z, [n, n]));
      
      title(sprintf('$r=%d, \\theta$=%0.2f', r, sh), 'Interpreter','latex')
      figSavePathName = sprintf('fourier_kernel r_%d shape_%dby100.png', r, 100*sh);
      saveas(hFig, figSavePathName)
    end
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

xun = mod(xun+0.3, 1);
K = prod(1 + (a)*constMult*bernPloy(xun),2);

end
