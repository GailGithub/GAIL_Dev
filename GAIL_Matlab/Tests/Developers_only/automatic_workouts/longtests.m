% LONGTESTS Drives all lengthy doctests, unit tests, workouts, and scripts

[GAILPATH,GAILVERSION,PATHNAMESEPARATOR,MATLABVERSION]  = GAILstart(0);
filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
    'gail_workouts-', datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.txt');
diary(filename)
tic;

format short

%% Workouts
% meanMC_g
doctest dt_meanMC_g_TrafficModel
run_handle('Test_meanMC_g')

% meanMCBer_g
 
% cubMC_g
format short
doctest dt_cubMC_g
run_handle('Test_cubMC_g')

% integral_g
run('ut_workout_integral_g')

% funappx_g
if MATLABVERSION >= 8
    run(ut_workoutfunappx_g)
    run(ut_convtest_funappx_g)
end

% funmin_g
if MATLABVERSION >= 8
    run(ut_workout_funmin_g)
end


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
%run_handle('DisplayTestResults_BlacknColor')

% meanMCBer_g paper
%run_handle('PlotmeanMCBernoulli_gResults')
run_handle('PlotRatioHoeffCLT')

% cubSobol_g paper
run_handle('RunTestCubatureonGeoAsianCallSobol');
run_handle('RunTestCubatureonKeisterSobol')

% cubLattice_g paper
run_handle('RunTestCubatureonGeoAsianCallLattice');
run_handle('RunTestCubatureonKeisterLattice');


% Function minimization thesis
try
    if usejava('jvm')
        run_handle('UniFunMin_Plot_Bump')
        run_handle('UniFunMin_Plot_Flat')
        run_handle('UniFunMin_Plot_TwoExtreme')
    end
    if MATLABVERSION >= 8
        run(ut_thesis_funmin01)
    end
catch
    display('Test for Papers/UniFunMin is wrongly coded. We skip it.')
    %fprintf(fid,'Test for Papers/UniFunMin is wrongly coded. We skip it.\n');
end

try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_meanMC_g);
    results=run(ut_meanMC_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Test ut_meanMC_g is wrongly coded. We skip it.')
    %fprintf(fid,'Test ut_meanMC_g is wrongly coded. We skip it.\n');
end

try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_meanMCBer_g);
    results=run(ut_meanMCBer_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Test ut_meanMCBer_g is wrongly coded. We skip it.')
    %fprintf(fid,'Test ut_meanMCBer_g is wrongly coded. We skip it.\n');
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
        % for i=1:size(failed,2)
        %    fprintf(fid,'%s\n',Tests(failed(i)).Name);
        % end
    end
catch
    display('Test ut_funappx01_g is wrongly coded. We skip it.')
end

warning('off','MATLAB:integral01_g:peaky')
try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_integral01_g);
    results=run(ut_integral01_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %    fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Test ut_integral01_g is wrongly coded. We skip it.')
end
warning('on','MATLAB:integral01_g:peaky')

try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_funmin01_g);
    results=run(ut_funmin01_g)
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %    fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Test ut_funmin01_g is wrongly coded. We skip it.')
end

%   try
%     Tests = matlab.unittest.TestSuite.fromClass(?ut_integralNoPenalty_g);
%     results=run(ut_integralNoPenalty_g)
%     if sum([results.Failed])>0
%       failed=find([results.Failed]>0);
%       for i=1:size(failed,2)
%         fprintf(fid,'%s\n',Tests(failed(i)).Name);
%       end
%     end
%   catch
%     display('Test ut_integralNoPenalty_g is wrongly coded. We skip it.')
%   end

time=toc;
disp(time)

diary off
