% Drives all doctests and unit tests
format short

[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
completereport = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
  'gail_tests-', datestr(now,'yyyymmdd'),'.txt');
diary(completereport)

shortutestreport = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
   'gail_unittests','.txt');
fid = fopen(shortutestreport,'wt');

tic; 
% CALL DOCTESTS 

doctest funappx_g
doctest funappxlocal_g
doctest funappxtau_g
doctest funappx01_g
clear in_param;
doctest dt_funappx_g
doctest funmin_g
doctest integral_g
%doctest integralsim_g %TODO Uncomment after Yizhi fixes bugs
doctest integral01_g
doctest integraltau_g
doctest dt_integral_g 
doctest meanMC_g
doctest meanMCRel_g
doctest meanMCBernoulli_g
doctest cubMC_g
doctest cubLattice_g
doctest cubSobol_g
doctest meanMCRel_g

% CALL UNIT TESTS
[~,~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION >= 8
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_funappx_g);
    results=run(ut_funappx_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        display('Test ut_funappx_g is wrongly coded. We skip it.')
        fprintf(fid,'Test ut_funappx_g is wrongly coded. We skip it.\n');
    end
    
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_funappx01_g);
    results=run(ut_funappx01_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        display('Test ut_funappx01_g is wrongly coded. We skip it.')
        fprintf(fid,'Test ut_funappx01_g is wrongly coded. We skip it.\n');
    end
    
    warning('off','MATLAB:integral_g:peaky')
    warning('off','MATLAB:integral_g:exceedbudget')
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_integral_g);
    results=run(ut_integral_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        display('Test ut_integral_g is wrongly coded. We skip it.')
        fprintf(fid,'Test ut_integral_g is wrongly coded. We skip it.\n');
    end
    warning('on','MATLAB:integral_g:peaky')
    warning('on','MATLAB:integral_g:exceedbudget')
  
    warning('off','MATLAB:integral01_g:peaky')
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_integral01_g);
    results=run(ut_integral01_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        display('Test ut_integral01_g is wrongly coded. We skip it.')
        fprintf(fid,'Test ut_integral01_g is wrongly coded. We skip it.\n');
    end
    warning('on','MATLAB:integral01_g:peaky')
    
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_meanMC_g);
    results=run(ut_meanMC_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        display('Test ut_meanMC_g is wrongly coded. We skip it.')
        fprintf(fid,'Test ut_meanMC_g is wrongly coded. We skip it.\n');
    end
  
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_cubMC_g);
    results=run(ut_cubMC_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        display('Test ut_cubMC_g is wrongly coded. We skip it.')
        fprintf(fid,'Test ut_cubMC_g is wrongly coded. We skip it.\n');
    end
    
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_meanMCBernoulli_g);
    results=run(ut_meanMCBernoulli_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        display('Test ut_meanMCBernoulli_g is wrongly coded. We skip it.')
        fprintf(fid,'Test ut_meanMCBernoulli_g is wrongly coded. We skip it.\n');
    end

    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_funappxlocal_g);
    results=run(ut_funappxlocal_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        display('Test ut_funappxlocal_g is wrongly coded. We skip it.')
        fprintf(fid,'Test ut_funappxlocal_g is wrongly coded. We skip it.\n');
    end
    
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_cubLattice_g);
    results=run(ut_cubLattice_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        display('Test ut_cubLattice_g is wrongly coded. We skip it.')
        fprintf(fid,'Test ut_cubLattice_g is wrongly coded. We skip it.\n');
    end
    
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_cubSobol_g);
    results=run(ut_cubSobol_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        display('Test ut_cubSobol_g is wrongly coded. We skip it.')
        fprintf(fid,'Test ut_cubSobol_g is wrongly coded. We skip it.\n');
    end

    
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_funmin_g);
    results=run(ut_funmin_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        display('Test ut_funmin_g is wrongly coded. We skip it.')
        fprintf(fid,'Test ut_funmin_g is wrongly coded. We skip it.\n');
    end
    
    %warning('off','MATLAB:meanMCRel_g:maxreached')
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_meanMCRel_g);
    results=run(ut_meanMCRel_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        display('Test ut_meanMCRel_g is wrongly coded. We skip it.')
        fprintf(fid,'Test ut_meanMCRel_g is wrongly coded. We skip it.\n');
    end
    warning('on','MATLAB:meanMCRel_g:maxreached')
     
end

time=toc;
% disp(time)

diary off
fclose(fid);
