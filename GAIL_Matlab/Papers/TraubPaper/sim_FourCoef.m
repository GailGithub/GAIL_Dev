% Simulations:
% - Complex Exponential
% - Direct Fourier sampling (no f'n evals)

%% Simulation settings
close all
clearvars
format short e
d = 3;     %dimension
coord_wt_vec = 2.^(0:-1:-d+1);
wt_pow = 8;
wv_num_max = 20;
wv_num_mtx = permn(-wv_num_max:wv_num_max,d); %compute this just once for all
nBasis = size(wv_num_mtx,1); %number of basis elements
n_app = 2^11; %number of points to use for approximating L_2 norm
% Assume that
%  x is nX x d
%  k is nK x d
u_wt_fun = @(k) prod(2.^((k~=0)/2) ./ max(1,coord_wt_vec .* abs(k)).^wt_pow ,2); %nK x 1
%v_wt_fun = @(k) -sign(k(:,1)) .* 2.^(sum(k~=0,2)/2); %nK x 1
v_wt_fun = @(k) -sign(k(:,1));
basisName = 'CompExp'; %complex Exponential
% inp_basis_fun = @(k,x) exp((2*pi*sqrt(-1))* (x * k'))./u_wt_fun(k)'; %nX x nK
% out_basis_fun = @(k,x) -exp((2*pi*sqrt(-1))* (x * k')); %nX x nK
inp_basis_fun = @(k,x) squeeze(prod(cos((2*pi)* (x .* reshape(k',1,d,size(k,1))) ...
   + reshape(k',1,d,size(k,1))*(pi/2)),2)) .* u_wt_fun(k)'; %nX x nK
out_basis_fun = @(k,x) squeeze(prod(cos((2*pi)* (x .* reshape(k',1,d,size(k,1))) ...
   + reshape(k',1,d,size(k,1))*(pi/2)),2)) .* v_wt_fun(k)'; %nX x nK
lambda = @(k) (2*pi)*(abs(k(:,1)).*u_wt_fun(k)); %nK x 1;
[~,wh_ordLambda] = sort(lambda(wv_num_mtx),'descend');

% Error tolerances
num_eps = 10; %no. of errors
min_log10_eps = -1;
max_log10_eps = 1;
eps_vec = 10.^(linspace(min_log10_eps,max_log10_eps,num_eps))';
eps_vec = flipud(eps_vec); %so that we visualize the smallest tolerance

%% Compute true function

% Compute Fourier coef and gammas for true function
%four_coef = -1 + 2*rand(nBasis,1); % uniform random numbers on [-1,1]
four_coef = randn(nBasis,1); % normal random numbers

p = sobolset(d);
p = scramble(p,'MatousekAffineOwen');
sob_pts = net(p,n_app); 
Sf_true = out_basis_fun(wv_num_mtx,sob_pts) * (four_coef .* lambda(wv_num_mtx));

%% Run algorithm for different error tolerances

err_vec(num_eps,1)=0; %container for errors
n_vec(num_eps,1)=0; %container for sample sizes
n_seq = [0 2.^(4:floor(log2(nBasis)))]';
n_n_seq = length(n_seq) - 1;
sig(n_n_seq,1) = 0;
a = 0.5;
inf_seq = a.^(0:1:n_n_seq);
inf_seq_nm = sqrt(1./(1 - a^2));

%Sum the Fourier coefficients
for kk = 1:n_n_seq
   range_wv_num = wh_ordLambda(n_seq(kk)+1:n_seq(kk+1));
   sig(kk) = norm(four_coef(range_wv_num) .* lambda(wv_num_mtx(range_wv_num,:)));
end
error_bd = sig*inf_seq_nm;
       

Sf_app(n_app,1) = 0;
for m = 1:length(eps_vec)
    %m
    
    %Compute partial sum
    
    % Algorithm:
    % 1) Compute sample size nn:
    wh_n = find(error_bd <= eps_vec(m),1,'first');
    n_vec(m) = n_seq(wh_n + 1);

    % 2) Compute true error between f and f_app    
    %evaluate f_app
    range_k = wh_ordLambda(1:n_vec(m));
    Sf_app = out_basis_fun(wv_num_mtx(range_k,:),sob_pts) * ...
       (four_coef(range_k) .* lambda(wv_num_mtx(range_k,:)));

    err_vec(m) = sqrt(mean((Sf_true - Sf_app).^2));

end

rat_vec = err_vec./eps_vec; %sample size vs error ratios
save(['sim_FourCoef_results_basis_' basisName '.mat']) %save results to plot later


[n_vec eps_vec rat_vec]


%%Plot results
%sim_FourCoef_plot_results(basisName)
 

