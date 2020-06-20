% plots the shift invariant fourier kernel for demo purposes
function plot_fourier_kernel()

n = 512;
shape_param = [0.2 0.8];
order = [1 2];
bvec = [0.2 0.8];
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
      Z(:,i) = kernel(xpts, 2*r, sh);
      leg_text{i} = sprintf('$r=%d,\\eta=%1.1f$', r, sh);
      i=i+1;
    end
  end

  plot(xpts,Z)
  xlabel('\(x\)')
  ylabel('\(C_{\theta}(x,0.3)\)')
  legend(leg_text, 'Interpreter','latex','location','best')
  axis tight

  figSavePathName = sprintf('fourier_kernel dim_1');
  gail.save_image(hFig, 'Paper_cubBayesLattice_g', figSavePathName)
end

if dim==2
  [X, Y] = meshgrid(l);
  
  for sh = shape_param
    for r = order
      hFig = figure('visible','on');
      set(hFig, 'units', 'inches', 'Position', [4 4 6.5 5.5])
      Z = kernel([X(:) Y(:)], 2*r, sh);
      meshc(X, Y, reshape(Z, [n, n]));
      
      title(sprintf('$r=%d, \\eta$=%0.2f', r, sh), 'Interpreter','latex')
      figSavePathName = sprintf('fourier_kernel r_%d shape_%dby100', r, 100*sh);
      gail.save_image(hFig, 'Paper_cubBayesLattice_g', figSavePathName)
    end
  end
  
  if true
    [X, Y] = meshgrid(l);
    
    for sh = shape_param
      for b = bvec
        hFig = figure('visible','on');
        set(hFig, 'units', 'inches', 'Position', [4 4 6.5 5.5])
        Z = kernelSmooth([X(:) Y(:)], b, sh);
        meshc(X, Y, reshape(Z, [n, n]));
        
        title(sprintf('$b=%0.2f, \\eta$=%0.2f', b, sh), 'Interpreter','latex')
        figSavePathName = sprintf('fourier_kernel b_%dby100 shape_%dby100', 100*b, 100*sh);
        %gail.save_image(hFig, 'Paper_cubBayesLattice_g', figSavePathName)
      end
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

function K = kernelSmooth(xun,b,shape)
  tempa = 2*pi*xun;  % Lattice
  tempb = b*cos(tempa);
  tempc = 1 + shape*(2*(tempb-1))./(1 + b^2 - 2*tempb);
  %tempc = 1 + shape*((b^2-1))./(1 + b^2 - 2*tempb);
  K = prod(tempc, 2);
end
