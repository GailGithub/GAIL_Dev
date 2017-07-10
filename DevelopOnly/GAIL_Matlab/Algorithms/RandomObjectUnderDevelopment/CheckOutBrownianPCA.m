% Try out Brownian PCA

format compact
clearvars

mauricioInp.timeDim.timeVector=0.1:0.1:1;
mauricioInp.payoffParam.optType={'amean'};
mauricioInp.priceParam.cubMethod = 'Sobol';
mauricioInp.priceParam.absTol = 0.001;

mauricioInp.bmParam.assembleType='PCA';
mauricio = optPrice(mauricioInp)

fredInp.timeDim.timeVector=0.1:0.1:1;
fredInp.payoffParam.optType={'amean'};
fredInp.priceParam.cubMethod = 'Sobol';
fredInp.priceParam.absTol = 0.001;
fred = optPrice(fredInp)

% tic, mauricioPaths = genPaths(mauricio,1e5); toc
% mean(mauricioPaths)
% var(mauricioPaths)
% 
% tic, fredPaths = genPaths(fred,1e5); toc
% mean(fredPaths)
% var(fredPaths)

tic, [price,out] = genOptPrice(mauricio), toc


tic, [price,out] = genOptPrice(fred), toc

