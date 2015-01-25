function GAIL_Publish()
% GAIL_PUBLISH  To generate html files in the GAIL subdirectory Documentation
if usejava('jvm')
    oldStatus = get(0,'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off')
    [GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION] = GAILstart(0);
    delete(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html',PATHNAMESEPARATOR,'*.png'))
    mfile_list = {'GAIL','funclist','help_funappx_g','help_funmin_g',...
        'help_integral_g', 'help_meanMC_g','help_meanMCBer_g', ...
        'help_cubMC_g','help_cubLattice_g','help_cubSobol_g','help_install'};
    for i=1:length(mfile_list),
        publish('help_install');
    end
    set(0, 'DefaultFigureVisible', oldStatus)
end

delete(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'gail_ug.m'))
delete(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html',PATHNAMESEPARATOR,'gail_ug.*'))
cat_cmd = 'cat ';
for i=1:length(mfile_list),
    cat_cmd = strcat([cat_cmd, ' ', GAILPATH,'Documentation',PATHNAMESEPARATOR,mfile_list{i},'.m', ' ']);
end
cat_cmd = strcat(cat_cmd, '>> gail_ug.m');
system(cat_cmd);
publish('gail_ug','pdf');
end