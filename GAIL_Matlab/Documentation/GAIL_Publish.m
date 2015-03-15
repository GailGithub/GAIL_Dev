function GAIL_Publish(ifGenerateHtml, ifGenerateLateX, ifBuildSearchIndex)
% GAIL_PUBLISH  To generate html files in the GAIL subdirectory Documentation
if usejava('jvm')

    oldStatus = get(0,'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off')
    [GAILPATH,GAILVERSION,PATHNAMESEPARATOR,~] = GAILstart(0);
    mfile_list = {'GAIL','help_license','help_readme','funclist',...
        'help_funappx_g','help_funmin_g',...
        'help_integral_g', 'help_meanMC_g','help_meanMCBer_g', ...
        'help_cubMC_g','help_cubLattice_g','help_cubSobol_g'};
    
    %% generate GAIL Documentation in HTML format
    if ifGenerateHtml
      delete(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html',PATHNAMESEPARATOR,'*.png'))
      for i=1:length(mfile_list),
        publish(mfile_list{i});
      end
    end
    
    %% generate GAIL Documentation in latex format
    if ifGenerateLateX
        s = computer;
        if all(s(1:2)=='PC') == 0
            delete(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html',PATHNAMESEPARATOR,'gail_ug.*'))
            cat_cmd = 'cat ';
            for i=1:length(mfile_list),
                cat_cmd = strcat([cat_cmd, ' ', GAILPATH,'Documentation',PATHNAMESEPARATOR,mfile_list{i},'.m', ' ']);
            end
            gailug_filename = strcat([GAILPATH,'Documentation',PATHNAMESEPARATOR,...
                'gail_ug',strrep(GAILVERSION, '.', '_'),'.m']);
            if exist(gailug_filename,'file') > 0
                delete(gailug_filename)
            end
            cat_cmd = strcat([cat_cmd, '>> ', gailug_filename]);
            system(cat_cmd);
            % publish(gailug_filename,'pdf');
            publish(gailug_filename,'latex');
        end
    end
    set(0, 'DefaultFigureVisible', oldStatus)
    warninfo = warning('query','MATLAB:doc:DocNotInstalled');
    warning('off', warninfo.identifier);
    
    %% build search index
    if ifBuildSearchIndex,
      builddocsearchdb(strcat(GAILPATH,'Documentation',PATHNAMESEPARATOR,'html'));
    end
    
    
    warning(warninfo.state,warninfo.identifier);
    fprintf('\nYou can go to help documentation ---> Supplemental Software to learn how to use GAIL.\n');
end
end
