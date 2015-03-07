function lattice_example
% In order to obtain the plotting properties in the article, we use
% Matlab's default properties that appear when starting it. If some
% properties have been changed, to obtain the same plots just restart
% Matlab again.

% set(0,'defaultLineMarkerSize',5) %large dots
[~,~,~,MATLABVERSION] = GAILstart(false); 
if usejava('jvm') || MATLABVERSION <= 7.12
a = 23;
b = 30;
%% Node set
h = [1 27 15];
n = 2^6;
j = (0:n-1)';
lattice = mod(j*h/n,1);
scatter3(lattice(:,1),lattice(:,2),lattice(:,3),'k','filled')
title('Node set')
view(a,b)
gail.save_eps('Paper_cubLattice_g', 'Paper_cubLattice_g_NodeSetExample');
figure
%% Dual Lattice
D = [];
for i = -10:10
    for j = -10:10
        for k = -10:10
            if sum(abs(mod(lattice*[i,j,k]',1)))==0
                D = [D ; [i,j,k]];
            end
        end
    end
end
scatter3(D(:,1),D(:,2),D(:,3),'k','filled')
title('Dual Lattice')
view(a,b)
gail.save_eps('Paper_cubLattice_g', 'Paper_cubLattice_g_DualSetExample');
close all
end