%%Plot results
function sim_FourCoef_plot_results(basisName)

load(['sim_FourCoef_results_basis_' basisName '.mat'], ...
   'rat_vec', 'eps_vec', 'n_vec', 'd', 'inp_basis_fun', 'out_basis_fun', ...
   'wv_num_mtx', 'wh_ordLambda', 'four_coef', 'lambda', ...
   'min_log10_eps', 'max_log10_eps') %load results for plotting

gail.InitializeDisplay %add some variables for nice plotting

figure
log10epsVec = log10(eps_vec);
h = scatter(rat_vec,n_vec,800,log10epsVec,'.'); %plot ratio of actual error to tolerance, with color corresonding to tolerance
hold on
set(gca,'XScale','log', 'YScale','log')
xlim([min([10.^floor(log10(rat_vec*0.8)); 0.1]) ...
 max([10.^(ceil(log10(rat_vec*1.2))); 1])])
ylim(10.^[floor(log10(min(n_vec))) ceil(log10(max(n_vec)))])
xlabel('\(||\partial^{\{1\}}f-\partial^{\{1\}}f_{\mbox{app}}||_{2}/\epsilon\)')
ylabel({'Sample size \(n\)'})
hcb = colorbar; %showing tolerance values
title(hcb,'\(\varepsilon\)','interpreter','latex')
tickVals = floor(min_log10_eps):ceil(max_log10_eps);
tickLabels = 10.^tickVals;
set(hcb,'Ticks',tickVals,'TickLabels',tickLabels, ...
   'Limits',[tickVals(1) tickVals(end)])
% [~,leg_icons] = legend(h, 'x^2', ...
%     'box','off','location','north','orientation','horizontal','interpreter','latex');
   
set(gcf,'Position',[200,200,1000,500]) %make figure big enough and the right aspect ratio
% leg_icons(2).Children.MarkerSize = 30; %make legend icons large enough
print -depsc SimFourCoefErr.eps

%   {'\(||\partial^{\{1\}}f-\partial^{\{1\}}f_{\mbox{app}}||_{2}/\epsilon\)', ...
%
%% Visualize (only a 2-d projection)
xcoord = 1; % coordinate with the biggest w
ycoord = 2; % coordinate with the next biggest w
nom_Val = 1/2;
n_grd = 100;
x_grd = 0:(2/n_grd):1;
n_grd_val = length(x_grd);
[xx,yy] = meshgrid(x_grd);
n_Vis = n_grd_val.^2;
x_Vis = nom_Val*ones(n_Vis,d);
x_Vis(:,xcoord) = xx(:);
x_Vis(:,ycoord) = yy(:);

%Input funcitons
f_true_Vis = inp_basis_fun(wv_num_mtx,x_Vis) * (four_coef .* lambda(wv_num_mtx));
m = length(eps_vec);
range_k = wh_ordLambda(1:n_vec(m));
f_app_Vis = inp_basis_fun(wv_num_mtx(range_k,:),x_Vis) * ...
   (four_coef(range_k) .* lambda(wv_num_mtx(range_k,:)));


%Solutions
Sf_true_Vis = out_basis_fun(wv_num_mtx,x_Vis) * (four_coef .* lambda(wv_num_mtx));
Sf_app_Vis = out_basis_fun(wv_num_mtx(range_k,:),x_Vis) * ...
   (four_coef(range_k) .* lambda(wv_num_mtx(range_k,:)));

    
figure
rotate3d on
Sf_true_Vis = reshape(Sf_true_Vis,[n_grd_val n_grd_val]);
surf(x_grd,x_grd,Sf_true_Vis); shading interp
xlabel(['\(x_{' int2str(xcoord) '}\)'])
ylabel(['\(x_{' int2str(ycoord) '}\)'])
zlabel(['\(\partial^{\{1\}}f(x_{' int2str(xcoord) '}, x_{' int2str(ycoord) ...
   '}, ' num2str(nom_Val) ', \ldots)\)'])
print -depsc SimFourCoefSol.eps


figure
rotate3d on
Sf_app_Vis = reshape(Sf_app_Vis,[n_grd_val n_grd_val]);
surf(x_grd,x_grd,Sf_app_Vis); shading interp
xlabel(['\(x_{' int2str(xcoord) '}\)'])
ylabel(['\(x_{' int2str(ycoord) '}\)'])
zlabel(['\(\partial^{\{1\}}f_{\mbox{app}}(x_{' int2str(xcoord) '}, x_{' int2str(ycoord) ...
   '}, ' num2str(nom_Val) ', \ldots)\)'])
print -depsc SimDirectSolAppx.eps

figure
rotate3d on
err_Vis = Sf_true_Vis - Sf_app_Vis;
surf(x_grd,x_grd,err_Vis); shading interp
xlabel(['\(x_{' int2str(xcoord) '}\)'])
ylabel(['\(x_{' int2str(ycoord) '}\)'])
zlabel(['\(\partial^{\{1\}}f_{\mbox{err}}(x_{' int2str(xcoord) '}, x_{' int2str(ycoord) ...
   '}, ' num2str(nom_Val) ', \ldots)\)'])
print -depsc SimDirectSolErr.eps


figure
rotate3d on
f_true_Vis = reshape(f_true_Vis,[n_grd_val n_grd_val]);
surf(x_grd,x_grd,f_true_Vis); shading interp
xlabel(['\(x_{' int2str(xcoord) '}\)'])
ylabel(['\(x_{' int2str(ycoord) '}\)'])
zlabel(['\(f(x_{' int2str(xcoord) '}, x_{' int2str(ycoord) ...
   '}, ' num2str(nom_Val) ', \ldots)\)'])
print -depsc SimFourCoefInpFun.eps


figure
rotate3d on
f_app_Vis = reshape(f_app_Vis,[n_grd_val n_grd_val]);
surf(x_grd,x_grd,f_app_Vis); shading interp
xlabel(['\(x_{' int2str(xcoord) '}\)'])
ylabel(['\(x_{' int2str(ycoord) '}\)'])
zlabel(['\(f_{\mbox{app}}(x_{' int2str(xcoord) '}, x_{' int2str(ycoord) ...
   '}, ' num2str(nom_Val) ', \ldots)\)'])
print -depsc SimDirectInpFunAppx.eps

figure
rotate3d on
err_Vis = f_true_Vis - f_app_Vis;
surf(x_grd,x_grd,err_Vis); shading interp
xlabel(['\(x_{' int2str(xcoord) '}\)'])
ylabel(['\(x_{' int2str(ycoord) '}\)'])
zlabel(['\(f_{\mbox{err}}(x_{' int2str(xcoord) '}, x_{' int2str(ycoord) ...
   '}, ' num2str(nom_Val) ', \ldots)\)'])
print -depsc SimDirectInpFunErr.eps

figure(1)


