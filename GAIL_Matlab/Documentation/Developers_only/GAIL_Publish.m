function GAIL_Publish()
% GAIL_PUBLISH  To generate html files in the GAIL subdirectory Documentation
if usejava('jvm')
    oldStatus = get(0,'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off')
    [GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION] = GAILstart(0);
    mfile_list = {'GAIL','funclist','help_funappx_g','help_funmin_g',...
        'help_integral_g', 'help_meanMC_g','help_meanMCBer_g', ...
        'help_cubMC_g','help_cubLattice_g','help_cubSobol_g','help_install'};
    
    %% generate GAIL User Guide in HTML format
    delete(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html',PATHNAMESEPARATOR,'*.png'))
    for i=1:length(mfile_list),
        publish(mfile_list{i});
    end
    
    %% generate GAIL User Guide in PDF format
    s = computer;
    if prod(s(1:2)=='PC') == 0
        delete(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html',PATHNAMESEPARATOR,'gail_ug.*'))
        cat_cmd = 'cat ';
        for i=1:length(mfile_list),
            cat_cmd = strcat([cat_cmd, ' ', GAILPATH,'Documentation',PATHNAMESEPARATOR,mfile_list{i},'.m', ' ']);
        end
        gailug_filename = strcat([GAILPATH,'Documentation',PATHNAMESEPARATOR,...
            'gail_ug',strrep(GAILVERSION, '.', '_'),'.m']);
        delete(gailug_filename)
        cat_cmd = strcat([cat_cmd, '>> ', gailug_filename]);
        system(cat_cmd);
        publish(gailug_filename,'pdf');
        publish(gailug_filename,'latex');
        set(0, 'DefaultFigureVisible', oldStatus)
    end
    
    builddocsearchdb(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html'));
    fprintf('\nYou can go to help documentation ---> Supplemental Software to learn how to use GAIL.\n');
   
end
end
