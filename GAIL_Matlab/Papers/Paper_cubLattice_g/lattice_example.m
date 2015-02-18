function lattice_example
a = 23;
b = 30;
%% Node set
h = [1 27 15];
n = 2^6;
j = (0:n-1)';
lattice = mod(j*h/n,1);
scatter3(lattice(:,1),lattice(:,2),lattice(:,3),'k','filled')
title('Node set')
gail.save_eps('Paper_cubLattice_g', 'Paper_cubLattice_g_NodeSetExample');
view(a,b)
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
gail.save_eps('Paper_cubLattice_g', 'Paper_cubLattice_g_DualSetExample');
view(a,b)
close all