function GAIL_Publish(ifGenerateHtml, ifGenerateLateX, ifBuildSearchIndex)
% GAIL_PUBLISH  To generate formatted documentation in the GAIL subdirectory Documentation
%
% ifGenerateHtml: true to generate HTML documentation
% ifGenerateLateX: true to generate tex source documentation files
% ifBuildSearchIndex: true to build search index in MATLAB
%
% Example 1:
% To generate search index for HTML documentation, run:
%    gail.GAIL_Publish(1,0,1)
%
% Example 2:
% To generate search index for HTML documentation and also LaTeX source, run:
%    gail.GAIL_Publish(1,1,1)
%
[~,~,~] = GAILstart;
if usejava('jvm')
    %oldStatus = get(0,'DefaultFigureVisible');
    %set(0, 'DefaultFigureVisible', 'off')
    [GAILPATH,GAILVERSION] = GAILstart(0);
    mfile_list = {'GAIL','help_license','help_readme','help_ReleaseNotes','funclist','demolist',...
        'help_funappx_g','help_funmin_g','help_integral_g',...
        'help_meanMC_g', 'help_meanMC_CLT','help_cubMC_g',...%'help_cubMC_CLT',...
        'help_cubLattice_g','help_cubSobol_g','help_cubBayesLattice_g', 'help_cubBayesNet_g',...
        'demo_funappx_g','demo_funappx_g1', 'demo_funappx_g2',...
        'demo_funmin_g','demo_funmin_g1', 'demo_funmin_g2',...
        'demo_integral_g','demo_integral_g1','demo_meanMC_g',...
        'count_success', 'demo_cubMC_g','demo_cubSobol_g','demo_cubLattice_g',...
        'demo_normal_probabilities','demo_cubBayesLattice_g','demo_cubBayesNet_g',...
        'demo_normal_probabilities_small','demo_meanMC_CLT'};
    wofile_list = {}; %{'Test_cubSobol_g'}; 
    
    %% generate GAIL Documentation in HTML format
    if ifGenerateHtml
      opts.stylesheet = strcat(GAILPATH,'Documentation',filesep,'mxdom2mathjaxbigfont.xsl');
      delete(strcat(GAILPATH,'Documentation',filesep,'html',filesep,'*.png'))
      for i=1:length(mfile_list)
        publish(mfile_list{i}, opts);
      end
      for i=1:length(wofile_list)
        publish(wofile_list{i}, opts);
        wopath = which(wofile_list{i});
        htmlpath = strrep(wopath, strcat([wofile_list{i},'.m']), strcat(['html',filesep]));
        %htmlfile = strcat([wofile_list{i}, '.html']);
        copyfile(fullfile(htmlpath,'*'),fullfile(GAILPATH,'Documentation','html'));
        %rmdir(htmlpath);
      end
    end
   
    
    %% build search index
    if ifBuildSearchIndex
        warninfo = warning('query','MATLAB:doc:DocNotInstalled');
        warning('off', warninfo.identifier);
        builddocsearchdb(strcat(GAILPATH,'Documentation',filesep,'html'));
        warning(warninfo.state,warninfo.identifier);
    end
    
    fprintf('\nYou can go to help documentation ---> Supplemental Software to learn how to use GAIL.\n');
    
    %% generate GAIL Documentation in latex format
    if ifGenerateLateX
        s = computer;
        if all(s(1:2)=='PC') == 0
            delete(strcat(GAILPATH,'Documentation',filesep,'Developers_only',filesep,'gail_ug.*'))
            cat_cmd = 'cat ';
            for i=1:length(mfile_list)
                cat_cmd = strcat([cat_cmd, ' ', GAILPATH,'Documentation',filesep,mfile_list{i},'.m', ' ']);
            end
            for i=1:length(wofile_list)
                wopath = which(wofile_list{i});
                wopath = strrep(wopath, strcat([wofile_list{i},'.m']), strcat(''));
                cat_cmd = strcat([cat_cmd, ' ', wopath, wofile_list{i},'.m' ]);
            end
            gailug_file = strcat(['gail_ug', strrep(GAILVERSION, '.', '_'),'.m']);
            gail_dir = strcat([GAILPATH, 'Documentation', filesep, 'Developers_only', filesep]);
            gailug_path = strcat([gail_dir, gailug_file]);
            if exist(gailug_path,'file') > 0
                delete(gailug_path)
            end
%           echo_cmd = strcat(['echo ', '"', 'function ', strrep(gailug_file,'.m',''), '"', ' > ',gailug_file]);
%           [status,result] = system(echo_cmd);
            cat_cmd = strcat([cat_cmd, '>> ', gailug_path]);
            [~,~] = system(cat_cmd);
%           echo_cmd = strcat(['echo ', '"', 'end', '"', ' >> ', gailug_file]);
%           [status,result] = system(echo_cmd);
%           publish(gailug_filename,'pdf');
            publish(gailug_path,'latex');
        end
    end
    %set(0, 'DefaultFigureVisible', oldStatus)
end
end
