%% Parameter settings
set(0,'defaultaxesfontsize',12,'defaulttextfontsize',12) %make font larger
set(0,'defaultLineLineWidth',1) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLegendInterpreter','latex') %latex axis labels
set(0,'defaultLineMarkerSize',7) %latex axis labels

m = 5;
n = 2^m;
%% Digital net
sobstr = sobolset(2);
sobstr = scramble(sobstr,'MatousekAffineOwen');
points = sobstr(1:n, :);

a = 0;
b = 0;
for i = 0:5;
    ind = 1:2^(i);
    scatter(points(ind,1),points(ind,2),'filled','k','marker','o')
    axis([0 1 0 1])
    grid on
    grid minor
    set(gca,'xminortick','on','xtick',[0:2^(-3):1])
    set(gca,'yminortick','on','ytick',[0:2^(-3):1])
    hold on
    % Plotting all grid lines for stratification
    for k = 1:a
        for r = 1:2^a-1;
            plot(r*[2^-a 2^-a], [0 1], 'r', 'LineWidth', 3)
        end
    end
    for k = 1:b
        for r = 1:2^b-1;
            plot([0 1], r*[2^-b 2^-b],  'r', 'LineWidth', 3)
        end
    end
    if a == b
        a = a + 1;
    else
        b = b + 1;
    end
    hold off
    if i == 5
        title(strcat('$\mathcal{P}_', num2str(i),'\oplus\mathbf{\Delta}$'), 'Interpreter', 'Latex')
        print('-depsc', strcat('sobol_p_', num2str(i),'.eps'))
    end
end

%% Rank-1 lattice
z1 = [1 31];
z2 = [1 7];

points1 = mod(bsxfun(@times, z1/n, [0:n-1]') + rand, 1);
points2 = mod(bsxfun(@times, z2/n, [0:n-1]') + rand, 1);
figure 
scatter(points1(:, 1), points1(:, 2), 'filled', 'k', 'marker', 'o')
axis([0 1 0 1])
title(strcat('$\mathcal{P}_', num2str(m),'\oplus\mathbf{\Delta}$'), 'Interpreter', 'Latex')
print('-depsc', strcat('lattice1_p_', num2str(m),'.eps'))
grid on
grid minor
figure 
scatter(points2(:, 1), points2(:, 2), 'filled', 'k', 'marker', 'o')
axis([0 1 0 1])
title(strcat('$\mathcal{P}_', num2str(m),'\oplus\mathbf{\Delta}$'), 'Interpreter', 'Latex')
print('-depsc', strcat('lattice2_p_', num2str(m),'.eps'))
grid on
grid minor
