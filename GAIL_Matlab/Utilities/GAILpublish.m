%GAILPUBLISH  To generate html files in the GAIL subdirectory Documentation
set(0, 'DefaultFigureVisible', 'off')
publish('GAIL');
publish('funclist');
publish('help_funappx_g');
publish('help_funappxlocal_g');
publish('help_integral_g');
publish('help_meanMCabs_g');
%publish('help_meanMCBernoulli_g');
publish('help_cubMCabs_g');
set(0, 'DefaultFigureVisible', 'on')
