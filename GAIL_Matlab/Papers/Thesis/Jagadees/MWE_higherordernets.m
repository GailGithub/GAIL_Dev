% Min working example to test higher order nets from Dirk buyens
function MWE_higherordernets()
format short
close all
set(0,'defaultaxesfontsize',18,'defaulttextfontsize',18, ... %make font larger
      'defaultLineLineWidth',3, ... %thick lines
      'defaultLineMarkerSize',25) %big dots
    
% fName = 'periodic'; f = @(x) prod(6*x.*(1-x),2); intf = 1; 
% fName = 'quad'; f = @(x) prod(2*x, 2); intf = 1; 
fName='Keister'; integrand = @(x,dim) keisterFunc(x,dim,1/sqrt(2)); exactIntf = @(dim) Keistertrue(dim);
ptransform = 'C1sin';

nvec = 2.^(1:18);

% Generators were downloaded from url:
% https://bitbucket.org/dnuyens/qmc-generators/src/3b059b39685a/DIGSEQ/?at=master
gen_names = {'exod2_base2_m20', 'nxmats/nx_s5_alpha3_m32', 'sobolmats/sobol_alpha3_Bs', };
figH = figure(); fc=1;
fig_size = [1 1 15 8];  % [1 1 9 6]
set(figH, 'units', 'inches', 'Position', fig_size)

for gen_name=gen_names(1:3)
  gen_name = gen_name{1};
  if strcmp(gen_name, 'exod2_base2_m20')
    genmat = gen_name;
  else
    genmat = load(['dirk_nuyens/qmc-generators/DIGSEQ/' gen_name '.col']);
  end
  for dim=[1 2 3]
    shift = rand(1,dim);
    f = @(x) integrand(x,dim);
    f = cubBayesLattice_g.doPeriodTx(f, ptransform) ;
    intf = exactIntf(dim);
    [cubf,errf] = compInteg(f,dim,nvec,genmat,intf,shift);
    outS.cubf = cubf; outS.errf = errf; outS.dim=dim; 
    % figure(); plot(xpts, f(xpts), '.')
    subplot(3,3,fc); fc=fc+1;
    loglog(nvec,errf,'.')
    axis(10.^[0 7 -16 0])
    set(gca,'Xtick',10.^[0 2 5 7], 'YTick',10.^[-16 -12 -8 -4 0])
    grid on
    g = strsplit(gen_name, '/'); title(sprintf('d=%d %s', dim, g{end}),'Interpreter','none')
  end
end
figSavePathName = sprintf('%s_%d_%s.png', fName, nvec(end), ptransform);
saveas(figH, figSavePathName)

fprintf('done')
end

% function [cubf,errf] = compInteg(f,dim,nvec,genmat,intf)
% nn = length(nvec);
% errf = zeros(nn,1);
% for i=1:nn
%   n=nvec(i);
%   digitalseq_b2g('init0', genmat)
%   xpts = digitalseq_b2g(dim, n)';
%   %plot(xpts,zeros(n,1),'.')
%   cubf = mean(f(xpts));
%   errf(i) = abs(intf - cubf);
% end
% end

% % using extensible pointset
% function [cubf,errf] = compInteg(f,dim,nvec,genmat,intf)
% nn = length(nvec);
% errf = zeros(nn,1);
% digitalseq_b2g('init0', genmat)
% n=nvec(1);
% xpts = digitalseq_b2g(dim, n)';
% for i=1:nn
%   %plot(xpts,zeros(n,1),'.')
%   cubf = mean(f(xpts));
%   errf(i) = abs(intf - cubf);
%   xpts = [xpts; digitalseq_b2g(dim, n)'];
%   n = 2*n;
% end
% end

% 
% % using extensible pointset
% function [cubf,errf] = compInteg(f,dim,nvec,genmat,intf)
% nn = length(nvec);
% errf = zeros(nn,1);
% digitalseq_b2g('init0', genmat)
% n=nvec(1);
% cubf=0;
% xpts = digitalseq_b2g(dim, n)';
% for i=1:nn
%   %plot(xpts,zeros(n,1),'.')
%   cubf = (cubf + mean(f(xpts)))/2;
%   errf(i) = abs(intf - cubf);
%   xpts = digitalseq_b2g(dim, n)';
%   n = 2*n;
% end
% end

% using extensible pointset
function [cubf,errf] = compInteg(f,dim,nvec,genmat,intf,shift)
nn = length(nvec);
errf = zeros(nn,1);
n=nvec(1);
if strcmp(genmat, 'exod2_base2_m20')
  xpts = cubBayesLattice_g.simple_lattice_gen(n,dim,shift,true);
else
  digitalseq_b2g('init0', genmat)
  xpts = digitalseq_b2g(dim, n)';
end
for i=1:nn
  %plot(xpts,zeros(n,1),'.')
  cubf = mean(f(xpts));
  errf(i) = abs(intf - cubf);
  if strcmp(genmat, 'exod2_base2_m20')
    n = 2*n;
    xpts_next = cubBayesLattice_g.simple_lattice_gen(n,dim,shift,false);
  else
    xpts_next = digitalseq_b2g(dim, n)'; n = 2*n;
  end
  xpts = [xpts; xpts_next];
end
end
