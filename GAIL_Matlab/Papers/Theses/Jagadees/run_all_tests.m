
!synclient HorizEdgeScroll=0 HorizTwoFingerScroll=0 %% disable horizontal scrolling
% gail.InitializeWorkspaceDisplay %initialize the workspace and the display parameters
format long

figSavePath = '/home/jagadees/MyWriteup/May4thweek_nonzero_mean/';
visiblePlot=false;

tstart=tic;
pdTx = {'C1sin', 'C2sin', 'Baker', 'C0', 'C1', 'none'};
for tx=pdTx
  for dim=[4 3 2]
    for bern=[2 4]
      TestExpCosBayesianCubature(dim,bern,tx{1},figSavePath,visiblePlot)
      if dim~=4
        TestMVN_BayesianCubature(dim,bern,tx{1},figSavePath,visiblePlot)
      end
      TestKeisterBayesianCubature(dim,bern,tx{1},figSavePath,visiblePlot)
    end
  end
end
toc(tstart)


if 0
cp MyWriteup/Apr1stweek/MVN/C2sin/* Dropbox/writeup/BeamerPresent/figures/ &&
cp MyWriteup/Apr1stweek/Keister/C2sin/* Dropbox/writeup/BeamerPresent/figures/ && 
cp MyWriteup/Apr1stweek/Exp\(cos\)/C1sin/* Dropbox/writeup/BeamerPresent/figures/

cp MyWriteup/Apr1stweek/Exp\(cos\)/C2sin/* Dropbox/writeup/BeamerPresent/figures/
end


error 'finished'
