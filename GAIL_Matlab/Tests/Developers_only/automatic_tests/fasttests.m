% Drives all doctests and unit tests
format short

GAILPATH = GAILstart(0);
completereport = strcat(GAILPATH,'OutputFiles',filesep,...
    'gail_tests-', datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.txt');
diary(completereport)

shortutestreport = strcat(GAILPATH,'OutputFiles',filesep,...
    'gail_unittests','.txt');
fid = fopen(shortutestreport,'wt');

tic

%% CALL DOCTESTs for individual algorithms moved to own files
fasttests_cubMC_g
fasttests_funappx_g
fasttests_integral_g

%% CALL DOCTESTS
tic; doctest gail.gail_in_param; time=toc
tic; doctest gail.gail1D_in_param; time=toc
tic; doctest gail.funmin_g_in_param; time=toc
tic; doctest gail.gailMD_in_param; time=toc

tic; doctest funmin_g; time=toc
tic; doctest meanMC_g; time=toc
tic; doctest cubLattice_g; time=toc
tic; doctest cubSobol_g; time=toc

tic; doctest assetPath; time=toc
tic; doctest optPayoff; time=toc
tic; doctest optPrice; time=toc
tic; doctest brownianMotion; time=toc
tic; doctest stochProcess; time=toc
tic; doctest whiteNoise; time=toc

%% CALL UNIT TESTS
[~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION < 8.1
    warning('Cannot run unit tests in MATLAB version before 8.1');
else
    
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
        display('Error: Test ut_funmin_g is wrongly coded. We skip it.')
        fprintf(fid,'Error: Test ut_funmin_g is wrongly coded. We skip it.\n');
    end
    
    try
        Tests = matlab.unittest.TestSuite.fromClass(?ut_funmin_g_end);
        results=run(ut_funmin_g_end)
        if sum([results.Failed])>0
            failed=find([results.Failed]>0);
            for i=1:size(failed,2)
                fprintf(fid,'%s\n',Tests(failed(i)).Name);
            end
        end
    catch
        display('Error: Test ut_funmin_g_end is wrongly coded. We skip it.')
        fprintf(fid,'Error: Test ut_funmin_g_end is wrongly coded. We skip it.\n');
    end
    
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
        display('Error: Test ut_meanMC_g is wrongly coded. We skip it.')
        fprintf(fid,'Error: Test ut_meanMC_g is wrongly coded. We skip it.\n');
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
        display('Error: Test ut_cubLattice_g is wrongly coded. We skip it.')
        fprintf(fid,'Error: Test ut_cubLattice_g is wrongly coded. We skip it.\n');
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
        display('Error: Test ut_cubSobol_g is wrongly coded. We skip it.')
        fprintf(fid,'Error: Test ut_cubSobol_g is wrongly coded. We skip it.\n');
    end
     
end

time=toc
% disp(time)

diary off
fclose(fid);
format
