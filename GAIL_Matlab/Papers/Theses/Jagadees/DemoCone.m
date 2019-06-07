
% gail.InitializeWorkspaceDisplay %initialize the workspace and the display parameters
format long
GAIL_path = GAILstart(0);
logSavePath=[GAIL_path,'Papers', filesep, 'Theses',filesep, 'Jagadees', filesep];

dim = 1;
b = 1;
n = 16;
A_rough = 1;
ftrue = @(x) exp(sum(cos(2*pi*x),2));
exactInteg = besseli(0,1)^dim;

shift = [0.1 0.5];

% genr = [1, 433461];
genr = [1 7];
% dual space
if dim==1
  genr_d = [32]';
else
  genr_d = [4 4]'; 
end

gen_lattice = @(n) mod(bsxfun(@times,(0:1/n:1-1/n)',genr(1:dim)),1); % unshifted in direct order
xpts_un = gen_lattice(n);
xpts = mod(bsxfun(@plus,xpts_un,shift(1:dim)),1);  % shifted in direct order

fnoise = @(x) (1 - exp(2*pi*sqrt(-1)*x*genr_d(1:dim)));
freal = @(x) (ftrue(x) + A_rough*real(fnoise(x)));
A_nice = -0.2;
fnice = @(x) (ftrue(x) + A_nice*real(fnoise(x)));

assert(all((ftrue(xpts_un)-real(freal(xpts_un)))==0))

if dim==1
  figH3 = figure; 
  xplot = [0:0.001:1]';

  plot(xplot, fnice(xplot))
  hold on  
  plot(xpts_un(:,1), ftrue(xpts_un(:,1)), '+', 'MarkerSize',6)
  plot(xplot, freal(xplot))
    
  xlabel('\( x \)','Interpreter','latex')
  ylabel('\( f(x) \)', 'Interpreter','latex')
  legend({'$f_{{nice}}$', '$f_{{true}}$', '$f_{{peaky}}$'},...
    'Interpreter','latex', 'location','best') % R2018a and earlier
  
  saveas(figH3, sprintf('cone_bayes_f_real.png'))
end

if dim==2
  [X,Y] = meshgrid(0:0.01:1,0:0.01:1);
  Z = freal([X(:) Y(:)]);
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

mu_nice = exactInteg + A_nice;
mu_real = exactInteg + A_rough;
mu_true = exactInteg;

gauss = @(x) (1/sqrt(2*pi*s2_MLE))*exp(-((x-m_MLE).^2)/(2*s2_MLE));

xa = 0:0.1:3;
figH4 = figure; 
plot(xa, gauss(xa))
hold on
plot(mu_real, 0, '.')
plot(mu_nice, 0, '.')
plot(exactInteg, 0, '.')

xlabel('$\mu$','Interpreter','latex')
ylabel('$\bf{P}\mu|_{(\boldmath{f}=\boldmath{y})}$', 'Interpreter','latex')
% [mu_real, 0], [mu_real 0.5];
xa_norm = @(t) (t - xa(1))/(xa(end)-xa(1));
axis tight

% annotation('textarrow',[xa_norm(mu_real) xa_norm(mu_real)],[0 0.25],'String','$\mu_{real}$')
% annotation('textarrow',[xa_norm(mu_nice) xa_norm(mu_nice)],[0 0.25],'String','$\mu_{nice}$')
% annotation('textarrow',[xa_norm(mu_true) xa_norm(mu_true)],[0 0.5],'String','$\mu_{true}$')

text(mu_real,0.15,'$\downarrow\mu_{peaky}$','Interpreter','latex')
text(mu_nice,0.2,'$\downarrow\mu_{nice}$','Interpreter','latex')
text(exactInteg,0.05,'$\downarrow\mu_{true}$','Interpreter','latex')

saveas(figH4, sprintf('cone_bayes_mu_pdf.png'))
  
fprintf('')

  