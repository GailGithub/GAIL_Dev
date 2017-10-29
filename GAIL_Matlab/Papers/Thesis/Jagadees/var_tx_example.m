

f = @(x)cos(x);

t = [0:0.001:1]; 

%% https://www.mathworks.com/matlabcentral/fileexchange/27991-tight-subplot-nh--nw--gap--marg-h--marg-w-
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
% [ha, pos] = tight_subplot(1,3,[.01 .03],[.1 .01],[.01 .01]) ;
hFig = figure(15);
set(hFig, 'units', 'inches', 'Position', [4 4 8.5 5.5])

[ha, pos] = tight_subplot(1,4, ...       % num_rows, num_colums
                          [.01 .06], ...  % gap_v, gap_h
                          [.18 .13], ...  % marg_vert_botm, marg_vert_top
                          [.05 .01]) ;    % marg_horz_left, marg_horz_right


axes(ha(1)); plot(t, t - sin(2*pi*t)/(2*pi)); 
xlabel('\(t\)'); title('\(x\)'); axis tight
ylim([0, 2]); xlim([0, 1])
%set(gca, 'XTickLabel', [0 0.5 0.1]);
yticks( [0 1 2]);
xticks( [0 0.5 1]);


%yticks([10^(-18) 10^(-15) 10^(-10) 10^(-5) 10^(-1)])  %(10.^[-18:5:-2])

    
axes(ha(2)); plot(t, f(t) );
xlabel('\(x\)'); title('\(f(x)\)'); axis tight
ylim([0, 2]); xlim([0, 1])
%set(gca, 'XTickLabel', [0 0.5 0.1]);
xticks( [ 0.5 1]);
yticks( [0 1 2]);
%set(gca, 'YTickLabel', []);

axes(ha(3)); plot(t, f(t - sin(2*pi*t)/(2*pi)).*(1-cos(2*pi*t)));
xlabel('\(t\)'); title('\(\tilde{f}(t)\)');  axis tight
ylim([0, 2]); xlim([0, 1])
xticks( [ 0.5 1]);
yticks( [0 1 2]);
%set(gca, 'YTickLabel', []);

tf1 = @(t) ( cos(t-(1/(2*pi)).*sin(2*pi*t)).*2*pi.*sin(2*pi*t) - ...
  sin(t-(1/(2*pi)).*sin(2*pi*t)).*(1-cos(2*pi*t)).^2);

axes(ha(4)); plot(t, tf1(t));
xlabel('\(t\)'); title('\(\tilde{f}^1(t)\)');  axis tight
ylim([-7, 7]); xlim([0, 1])
xticks( [ 0.5 1]);
%set(gca, 'YTickLabel', []);

saveas(hFig, sprintf('var_tx_example.eps' ))

fprintf('done');


