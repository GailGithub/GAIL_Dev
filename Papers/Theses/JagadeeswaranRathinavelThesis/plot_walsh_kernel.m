
function plot_walsh_kernel()
n = 512;
shape_param = [0.2 0.5 0.8];
order = [1];
dim = 1;

kernelFunc = cubBayesNet_g.BuildKernelFunc(order);
      
% using uniform points
xpts = linspace(0,1,n)';
xpts = [0:1/n:1-1/n]';

if dim==1
  Z = zeros(size(xpts, 1), 3);
  hFig = figure('visible','on');
  set(hFig, 'units', 'inches', 'Position', [4 4 6.5 5.5])
  hold on
  i=1;
  leg_text = cell(3,1);
  for r = order
    for eta = shape_param
      Z(:,i) = prod(1 + eta*kernelFunc(mod(xpts, 1)),2);

      leg_text{i} = sprintf('$\\eta=%1.1f$', eta);
      i=i+1;
    end
  end

  plot(xpts,Z, '.', 'MarkerSize', 20)
  % line_styles = {'--',':','-.','-'};
  % pointShapes = {'p', 'o','h','>','d','^','v','p','h'};

  % for i=1:size(Z,2)
    % plot(xpts,Z(:,i), pointShapes{i}, 'MarkerSize', 3)
  % end
  
  xlabel('\(x\)')
  ylabel('\(C_{\theta}(x,0)\)')
  legend(leg_text, 'Interpreter','latex','location','best')
  axis tight

  figSavePathName = sprintf('walsh_kernel dim_1');
  gail.save_image(hFig, 'JagadeeswaranRathinavelThesis', figSavePathName)
  F = getframe(hFig);
  % imwrite(rgb2gray(F.cdata), 'walsh_kernel_dim_1-gray.png', 'png')  
end

end
