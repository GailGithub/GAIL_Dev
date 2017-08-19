function GAIL_Publish(ifGenerateHtml, ifGenerateLateX, ifBuildSearchIndex)
% GAIL_PUBLISH  To generate formatted documentation in the GAIL subdirectory Documentation
%
% ifGenerateHtml: true to generate HTML documentation
% ifGenerateLateX: true to generate tex source documentation files
% ifBuildSearchIndex: true to build search index in MATLAB
%
[GAILPATH,GAILVERSION,MATLABVERSION] = GAILstart;
if usejava('jvm')
    %oldStatus = get(0,'DefaultFigureVisible');
    %set(0, 'DefaultFigureVisible', 'off')
    [GAILPATH,GAILVERSION] = GAILstart(0);
    mfile_list = {'GAIL','help_license','help_readme','funclist','demolist',...
        'help_funappx_g','help_funmin_g',...
        'help_integral_g', 'help_meanMC_g', ...
        'help_cubMC_g','help_cubLattice_g','help_cubSobol_g',... % demos below
        'demo_funappx_g','demo_funappx_g1', 'demo_funappx_g2',...
        'demo_funmin_g','demo_funmin_g1', 'demo_funmin_g2',...
        'demo_integral_g','demo_integral_g1',...
        'demo_meanMC_g','count_success',...
        'demo_cubMC_g','demo_normal_probabilities_cubMC',...
        'demo_cubSobol_g','demo_normal_probabilities'};
    wofile_list = {'Test_cubSobol_g'}; 
    %% generate GAIL Documentation in HTML format
    if ifGenerateHtml
      opts.stylesheet = strcat(GAILPATH,'Documentation',filesep,'mxdom2mathjaxbigfont.xsl');
      delete(strcat(GAILPATH,'Documentation',filesep,'html',filesep,'*.png'))
      for i=1:length(mfile_list),
        publish(mfile_list{i}, opts);
      end
      for i=1:length(wofile_list),
        publish(wofile_list{i}, opts);
        wopath = which(wofile_list{i});
        htmlpath = strrep(wopath, strcat([wofile_list{i},'.m']), strcat(['html',filesep]));
        %htmlfile = strcat([wofile_list{i}, '.html']);
        copyfile(fullfile(htmlpath,'*'),fullfile(GAILPATH,'Documentation','html'));
        %rmdir(htmlpath);
      end
    end
    
    %% generate GAIL Documentation in latex format
    if ifGenerateLateX
        s = computer;
        if all(s(1:2)=='PC') == 0
            delete(strcat(GAILPATH,'Documentation',filesep,'html',filesep,'gail_ug.*'))
            cat_cmd = 'cat ';
            for i=1:length(mfile_list),
                cat_cmd = strcat([cat_cmd, ' ', GAILPATH,'Documentation',filesep,mfile_list{i},'.m', ' ']);
            end
            for i=1:length(wofile_list),
                wopath = which(wofile_list{i});
                wopath = strrep(wopath, strcat([wofile_list{i},'.m']), strcat(['']));
                cat_cmd = strcat([cat_cmd, ' ', wopath, wofile_list{i},'.m' ]);
            end
            gailug_file = strcat(['gail_ug',strrep(GAILVERSION, '.', '_'),'.m']);
            gailug_path = strcat([GAILPATH,'Documentation',filesep,gailug_file]);
            if exist(gailug_path,'file') > 0
                delete(gailug_path)
            end
%           echo_cmd = strcat(['echo ', '"', 'function ', strrep(gailug_file,'.m',''), '"', ' > ',gailug_file]);
%           [status,result] = system(echo_cmd);
            cat_cmd = strcat([cat_cmd, '>> ', gailug_file]);
            [status,result] = system(cat_cmd);
%           echo_cmd = strcat(['echo ', '"', 'end', '"', ' >> ', gailug_file]);
%           [status,result] = system(echo_cmd);
%           publish(gailug_filename,'pdf');
            publish(gailug_path,'latex');
        end
    end
    %set(0, 'DefaultFigureVisible', oldStatus)
    
    %% build search index
    if ifBuildSearchIndex,
        warninfo = warning('query','MATLAB:doc:DocNotInstalled');
        warning('off', warninfo.identifier);
        builddocsearchdb(strcat(GAILPATH,'Documentation',filesep,'html'));
        warning(warninfo.state,warninfo.identifier);
    end
    
    fprintf('\nYou can go to help documentation ---> Supplemental Software to learn how to use GAIL.\n');
end
end
