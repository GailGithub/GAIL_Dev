function DemoCone()
% gail.InitializeWorkspaceDisplay %initialize the workspace and the display parameters
% format long
GAIL_path = GAILstart(0);
logSavePath=[GAIL_path,'Papers', filesep, 'Theses',filesep, 'Jagadees', filesep];

dim = 1;
b = 1;
n = 4;
A_peak = 2.5;
ftrue = @(x) exp(sum(cos(2*pi*x),2));
exactInteg = besseli(0,1)^dim;

shift = [0.1 0.5];

% genr = [1, 433461];
genr = [1 7];
% dual space
if dim==1
  genr_d = [n]';
else
  genr_d = [4 4]';
end

gen_lattice = @(n) mod(bsxfun(@times,(0:1/n:1-1/n)',genr(1:dim)),1); % unshifted in direct order
xpts_un = gen_lattice(n);
xpts = mod(bsxfun(@plus,xpts_un,shift(1:dim)),1);  % shifted in direct order

fnoise = @(x) (1 - exp(2*pi*sqrt(-1)*x*genr_d(1:dim)));
fpeak = @(x) (ftrue(x) + A_peak*real(fnoise(x)));
% A_nice = -0.2;
% fnice = @(x) (ftrue(x) + A_nice*real(fnoise(x)));

assert(all((ftrue(xpts_un)-real(fpeak(xpts_un)))==0))

if dim==1
  figH3 = figure;
  set(figH3, 'units', 'inches', 'Position', [7 2 8 6])
  xplot = [0:0.001:1]';
  
  shape = 0.9;
  xun = xpts_un; %(1:2:end);
  [K,kvec] = kernelFun(xun,xun,shape);
  M = kernelFun(xplot,xun,shape);
  Kinvy = (K\ftrue(xun(:,1)));
  ynice = M*Kinvy;
  muhat = kvec'*Kinvy;
  
  % plot(xplot, fnice(xplot))
  plot(xplot, ynice, 'b')
  hold on
  plot(xplot, fpeak(xplot), 'r')
  plot(xpts_un(:,1), ftrue(xpts_un(:,1)), 'g+', 'MarkerSize',6)
  
  xlabel('\( x \)','Interpreter','latex')
  ylabel('\( f(x) \)', 'Interpreter','latex')
  legend({'$f_{\textup{nice}}$', '$f_{\textup{peaky}}$', ...
    '$f_{\textup{smooth}} = f_{\textup{nice}} = f_{\textup{peaky}}$'},...
    'Interpreter','latex', 'location','best') % R2018a and earlier
  
  gail.save_image(figH3, 'JagadeeswaranRathinavelThesis', sprintf('cone_bayes_f_real'))
end

if dim==2
  [X,Y] = meshgrid(0:0.01:1,0:0.01:1);
  Z = fpeak([X(:) Y(:)]);
  figure; surf(X,Y,reshape(Z, size(X)))
end

r = 2;
theta = 0.5;
constMult = -(-1)^(r/2)*((2*pi)^r)/factorial(r);
bernPoly = @(x)(-x.*(1-x) + 1/6);
kernel = @(x,theta) prod(1 + theta*constMult*bernPoly(x), 2);

lambda = real(fft(kernel(xpts_un,theta)));

y = ftrue(xpts_un);
ftilde = fft(y);

m_MLE = ftilde(1)/n;
s2_MLE = sum((abs(ftilde(2:end)).^2)./lambda(2:end))/(n^2);

mu_nice = muhat;  %exactInteg + A_nice;
mu_peak = exactInteg + A_peak;
mu_true = exactInteg;

gauss = @(x) (1/sqrt(2*pi*s2_MLE))*exp(-((x-m_MLE).^2)/(2*s2_MLE));

% xa = 0:0.001:3;
xa = -4:0.001:6;

figH4 = figure;
set(figH4, 'units', 'inches', 'Position', [7 2 8 5])

plot(xa, gauss((xa*sqrt(s2_MLE)+m_MLE)), 'k')
hold on
plot( (mu_peak-m_MLE)/sqrt(s2_MLE), 0, 'r.')
plot( (mu_nice-m_MLE)/sqrt(s2_MLE), 0, 'b.')
% plot( 0, 0, 'g.')

xlabel('$\sigma$','Interpreter','latex')
% ylabel('$\rho_\mu|_{(\boldmath{f}=\boldmath{y})}$', 'Interpreter','latex')
% [mu_peak, 0], [mu_peak 0.5];
% xa_norm = @(t) (t - xa(1))/(xa(end)-xa(1));
axis tight

% annotation('textarrow',[xa_norm(mu_peak) xa_norm(mu_peak)],[0 0.25],'String','$\mu_{real}$')
% annotation('textarrow',[xa_norm(mu_nice) xa_norm(mu_nice)],[0 0.25],'String','$\mu_{nice}$')
% annotation('textarrow',[xa_norm(mu_true) xa_norm(mu_true)],[0 0.5],'String','$\mu_{true}$')
legend({'$\rho_\mu|_{(\boldmath{f}=\boldmath{y})}$', ...
  '$\mu_\textup{peaky}$', '$\mu_\textup{nice}$'}, ...  %, '$\hat{\mu}$'
  'location','best','Interpreter','latex')
% temp = [{'${\vert\mu-\widehat{\mu} \vert}/{\varepsilon}=1$'}, temp'];
% legend(temp,'location','best','Interpreter','latex'); axis tight

%text(mu_peak,0.15,'$\downarrow\mu_{peaky}$','Interpreter','latex')
%text(mu_nice,0.2,'$\downarrow\mu_{nice}$','Interpreter','latex')
%text(exactInteg,0.05,'$\downarrow\mu_{true}$','Interpreter','latex')

gail.save_image(figH4, 'JagadeeswaranRathinavelThesis', sprintf('cone_bayes_mu_pdf'))

fprintf('')
end


function [K,kvec,k0] = kernelFun(x,t,shape)

[nx,d] = size(x);
[nt,~] = size(t);
domain = [zeros(1,d); ones(1,d)];

K = ones(nx,nt);

diffdom = diff(domain,1,1);
shdiffdom = shape*diffdom;
k0 = prod((- 6 + 4*shdiffdom +exp(-shdiffdom).*(6 + 2*shdiffdom))./shdiffdom.^2);
kvec = ones(nx,1)*(2^d/prod(shdiffdom));
for k = 1:d
  tempa = shape*abs(bsxfun(@minus,x(:,k),t(:,k)'));
  K = K.*exp(-tempa).*(1 + tempa);
  tempb = shape*(x(:,k)-domain(1,k));
  tempc = shape*(domain(2,k) - x(:,k));
  kvec = kvec.*(2 - exp(-tempc).*(1+tempc/2) ...
    - exp(-tempb).*(1+tempb/2));
end

end
