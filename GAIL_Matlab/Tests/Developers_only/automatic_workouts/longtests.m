% LONGTESTS Drives all lengthy doctests, unit tests, workouts, and scripts

[GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION]  = GAILstart(0);
filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
  'gail_workouts-', datestr(now,'yyyymmdd'),'.txt');
diary(filename)
tic; 

format short
  
%% Workouts 
% meanMC_g
doctest dt_meanMC_g_TrafficModel
run_handle('Test_meanMC_g')

% meanMCBer_g
run_handle('Test_meanMCBer_g')

% cubMC_g
format short
doctest dt_cubMC_g
run_handle('Test_cubMC_g')

% integral_g
run('ut_workout_integral_g')

% funappx_g
if MATLABVERSION >= 8  
    run(ut_workoutfunappx_g)
end

% funmin_g
run_handle('workout_ErrXToleranceTest')	
run_handle('workout_XToleranceTest.m')
run_handle('workout_ErrToleranceTest.m')
run_handle('workout_TwoExtremeTest.m')


%% Papers
% Cone paper
if usejava('jvm')
  run_handle('ConesPaperFoolFunctions')
end
if MATLABVERSION >= 8  
    run(ut_ConesPaper)
end

% MCQMC paper
run_handle('RunTestcubMConGeoAsianCall')
run_handle('RunTestcubMConGaussian')
run_handle('RunTestcubMConGaussiand1')

% meanMCBer_g paper
if MATLABVERSION >= 8  
    run(ut_meanMCBer_g)
end

% cubSobol_g paper
run_handle('RunTestCubatureonGeoAsianCall')
run_handle('RunTestCubatureonKeister')

% cubLattice_g paper
run_handle('TestFourierTransform');
run_handle('test_cubLattice');
run_handle('RunTestCubatureonGeoAsianCall');
run_handle('RunTestCubatureonKeister')

% Function minimization
if usejava('jvm')
  run_handle('UniFunMin_Plot_Bump')
  run_handle('UniFunMin_Plot_Flat')
  run_handle('UniFunMin_Plot_TwoExtreme')
end
if MATLABVERSION >= 8  
    run(ut_UniFunMin)
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
    display('Test ut_meanMC_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_meanMC_g is wrongly coded. We skip it.\n');
  end

 try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_meanMCBer_g);
    results=run(ut_meanMCBer_g)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display('Test ut_meanMCBer_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_meanMCBer_g is wrongly coded. We skip it.\n');
 end
  
%% doctests and unit tests for deprecated algos
doctest funappxtau_g
doctest funappxglobal_g
doctest dt_funappxglobal_g
doctest funappx01_g
doctest funmin01_g
doctest integral01_g
doctest integraltau_g

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
    Tests = matlab.unittest.TestSuite.fromClass(?ut_funmin01_g);
    results=run(ut_funmin01_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
catch
    display('Test ut_funmin01_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_funmin01_g is wrongly coded. We skip it.\n');
end

time=toc;
disp(time)

diary off