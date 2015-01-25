function GAIL_Publish()
% GAIL_PUBLISH  To generate html files in the GAIL subdirectory Documentation
if usejava('jvm')
    oldStatus = get(0,'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off')
    [GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION] = GAILstart(0);
    delete(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html',PATHNAMESEPARATOR,'*.png'))
    publish('GAIL');
    publish('funclist');
    publish('help_funappx_g');
    publish('help_integral_g');
    publish('help_meanMC_g');
    publish('help_meanMCBer_g');
    publish('help_cubMC_g');
    publish('help_funmin_g');
    %publish('help_cubLattice_g');
    %publish('help_cubSobol_g');
    publish('help_install');
    set(0, 'DefaultFigureVisible', oldStatus)
end

delete(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'gail_ug.m'))
delete(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html',PATHNAMESEPARATOR,'gail_ug.*'))
system('cat GAIL.m funclist.m help_meanMC_g.m help_funappx_g.m help_integral_g.m GAIL_Publish.m help_cubMC_g.m help_funmin_g.m help_meanMCBer_g.m >> gail_ug.m');

publish('gail_ug','pdf');
end