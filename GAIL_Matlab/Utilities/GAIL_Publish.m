function GAIL_Publish()
% GAIL_PUBLISH  To generate html files in the GAIL subdirectory Documentation
set(0, 'DefaultFigureVisible', 'off')
publish('GAIL');
publish('funclist');
publish('help_funappx_g');
publish('help_integral_g');
publish('help_meanMC_g');
publish('help_meanMCBer_g');
publish('help_cubMC_g');
publish('help_funmin_g');
publish('help_cubLattice_g');
publish('help_cubSobol_g');
publish('help_install');
set(0, 'DefaultFigureVisible', 'on')
end