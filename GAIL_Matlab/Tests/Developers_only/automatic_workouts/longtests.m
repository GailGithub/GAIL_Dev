% LONGTESTS Drives all lengthy doctests, unit tests, workouts, and scripts

[GAILPATH,GAILVERSION,MATLABVERSION]  = GAILstart(false);
filename = strcat(GAILPATH,'OutputFiles',filesep,...
    'gail_workouts-', datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.txt');
diary(filename)
tic;

format short

%% utilties
run_handle_ut('ut_save_mat')
if usejava('jvm') || MATLABVERSION <= 7.12
  run_handle_ut('ut_save_eps')
end

%% Workouts
longtests_cubMC_g
longtests_funappx_g
longtests_integral_g

% funappx_g
%doctest par_funappx_g

% funmin_g
%doctest par_funmin_g

% meanMC_g
doctest dt_meanMC_g_TrafficModel

if 0 % moved out to logtests_cubMC_g
% cubMC_g
format short
doctest dt_cubMC_g
end

if 0
% integral_g
try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_workout_integral_g);
    results=run(ut_workout_integral_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Error: Test ut_workout_integral_g is wrongly coded. We skip it.')
    %fprintf(fid,'Error: Test ut_workout_integral_g is wrongly coded. We skip it.\n');
end
end

if 0 % moved out to logtests_funappx_g
% funappx_g
try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_workout_funappx_g);
    results=run(ut_workout_funappx_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Error: Test ut_workout_funappx_g is wrongly coded. We skip it.')
    %fprintf(fid,'Error: Test ut_workout_funappx_g is wrongly coded. We skip it.\n');
end
try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_convtest_funappx_g);
    results=run(ut_convtest_funappx_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Error: Test ut_convtest_funappx_g is wrongly coded. We skip it.')
    %fprintf(fid,'Error: Test ut_convtest_funappx_g is wrongly coded. We skip it.\n');
end

% try
%     Tests = matlab.unittest.TestSuite.fromClass(?ut_par_funappx_g);
%     results = run(ut_par_funappx_g);
%     if sum([results.Failed])>0
%         failed=find([results.Failed]>0);
%         for i=1:size(failed,2)
%             fprintf(fid,'%s\n',Tests(failed(i)).Name);
%         end
%     end
% catch
%     display('Error: Test ut_par_funappx_g is wrongly coded. We skip it.')
% end
end


% funmin_g
% try
%     Tests = matlab.unittest.TestSuite.fromClass(?ut_workout_funmin_g);
%     results=run(ut_workout_funmin_g);
%     if sum([results.Failed])>0
%         failed=find([results.Failed]>0);
%         %for i=1:size(failed,2)
%         %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
%         %end
%     end
% catch
%     display('Error: Test ut_workout_funmin_g is wrongly coded. We skip it.')
%     %fprintf(fid,'Error: Test ut_workout_funmin_g is wrongly coded. We skip it.\n');
% end

% try
%     Tests = matlab.unittest.TestSuite.fromClass(?ut_par_funmin_g);
%     results = run(ut_par_funmin_g);
%     if sum([results.Failed])>0
%         failed=find([results.Failed]>0);
%         for i=1:size(failed,2)
%             fprintf(fid,'%s\n',Tests(failed(i)).Name);
%         end
%     end
% catch
%     display('Error: Test ut_par_funmin_g is wrongly coded. We skip it.')
% end

% try
%     Tests = matlab.unittest.TestSuite.fromClass(?ut_par_funmin_g_end);
%     results = run(ut_par_funmin_g_end);
%     if sum([results.Failed])>0
%         failed=find([results.Failed]>0);
%         for i=1:size(failed,2)
%             fprintf(fid,'%s\n',Tests(failed(i)).Name);
%         end
%     end
% catch
%     display('Error: Test ut_par_funmin_g_end is wrongly coded. We skip it.')
% end

%cubQMC
format short
doctest dt_cubSobol_g
doctest dt_cubLattice_g

%% Papers
% Cone paper
if usejava('jvm') || MATLABVERSION <= 7.12
    run_handle('ConesPaperFoolFunctions')
end
try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_ConesPaper);
    results=run(ut_ConesPaper);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Error: Test ut_ConesPaper is wrongly coded. We skip it.')
    %fprintf(fid,'Error: Test ut_ConesPaper is wrongly coded. We skip it.\n');
end

if 0 % moved out to logtests_cubMC_g
% MCQMC paper
run_handle('MCQMC2012Figs')
run_handle('FoolAutomaticAlgorithms')
run_handle('RunTestcubMConGeoAsianCall')
run_handle('RunTestcubMConGaussian')
run_handle('RunTestcubMConGaussiand1')
try
    DisplayTestResults_BlacknColor({'ex1', 'ex2', 'ex3'},'black')
catch
    display('Error: DisplayTestResults_BlacknColor is wrongly coded. We skip it.')
    %fprintf(fid,'Error: DisplayTestResults_BlacknColor is wrongly coded. We skip it.\n');
end
end

% cubSobol_g paper
try
  SobolWalshPict;
catch
    display('Error: SobolWalshPict is wrongly coded. We skip it.')
end
try
  WalshFourierCoeffDecay;
catch
    display('Error: WalshFourierCoeffDecay is wrongly coded. We skip it.')
end
%run_handle('RunTestCubatureonGeoAsianCallSobol');
%run_handle('RunTestCubatureonKeisterSobol')
try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_Papers_cubSobol_g);
    results=run(ut_Papers_cubSobol_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Error: Test ut_Papers_cubSobol_g is wrongly coded. We skip it.')
    %fprintf(fid,'Error: Test ut_Papers_cubSobol_g is wrongly coded. We skip it.\n');
end

% cubLattice_g paper
try
  lattice_example;
catch
    display('Error: lattice_example is wrongly coded. We skip it.')
end
try
  FourierCoeffDecayPict;
catch
    display('Error: FourierCoeffDecayPict is wrongly coded. We skip it.')
end
%run_handle('RunTestCubatureonGeoAsianCallLattice');
%run_handle('RunTestCubatureonKeisterLattice');
try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_Papers_cubLattice_g);
    results=run(ut_Papers_cubLattice_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Error: Test ut_Papers_cubLattice_g is wrongly coded. We skip it.')
    %fprintf(fid,'Error: Test ut_Papers_cubLattice_g is wrongly coded. We skip it.\n');
end

% Function minimization thesis
% try
%     if usejava('jvm') || MATLABVERSION <= 7.12
%         run_handle('UniFunMin_Plot_Bump')
%         run_handle('UniFunMin_Plot_Flat')
%         run_handle('UniFunMin_Plot_TwoExtreme')
%     end
% catch
%     display('Error: Test for Papers/UniFunMin is wrongly coded. We skip it.')
%     %fprintf(fid,'Error: Test for Papers/UniFunMin is wrongly coded. We skip it.\n');
% end
% try
%     Tests = matlab.unittest.TestSuite.fromClass(?ut_thesis_funmin01_g);
%     results=run(ut_thesis_funmin01_g);
%     if sum([results.Failed])>0
%         failed=find([results.Failed]>0);
%         % for i=1:size(failed,2)
%         %    fprintf(fid,'%s\n',Tests(failed(i)).Name);
%         % end
%     end
% catch
%     display('Error: Test ut_thesis_funmin01_g is wrongly coded. We skip it.')
% end

%% doctests and unit tests for deprecated algos
doctest funappxtau_g
doctest funappxglobal_g
%doctest dt_funappxglobal_g
% doctest funmin01_g
doctest integral01_g
doctest integraltau_g
doctest meanMCabs_g
doctest cubMCabs_g;
doctest cubLattice_old_g;
doctest cubSobol_old_g;

% try
%     Tests = matlab.unittest.TestSuite.fromClass(?ut_funappxglobal_g);
%     results=run(ut_funappxglobal_g);
%     if sum([results.Failed])>0
%         failed=find([results.Failed]>0);
%         for i=1:size(failed,2)
%             fprintf(fid,'%s\n',Tests(failed(i)).Name);
%         end
%     end
% catch
%     display('Error: Test ut_funappxglobal_g is wrongly coded. We skip it.')
%     %fprintf(fid,'Error: Test ut_funappxglobal_g is wrongly coded. We skip it.\n');
% end

warning('off','GAIL:integral01_g:peaky')
try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_integral01_g);
    results=run(ut_integral01_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %    fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Error: Test ut_integral01_g is wrongly coded. We skip it.')
end
warning('on','GAIL:integral01_g:peaky')

% try
%     Tests = matlab.unittest.TestSuite.fromClass(?ut_funmin01_g);
%     results=run(ut_funmin01_g);
%     if sum([results.Failed])>0
%         failed=find([results.Failed]>0);
%         %for i=1:size(failed,2)
%         %    fprintf(fid,'%s\n',Tests(failed(i)).Name);
%         %end
%     end
% catch
%     display('Error: Test ut_funmin01_g is wrongly coded. We skip it.')
% end

try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_meanMCabs_g);
    results=run(ut_meanMCabs_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
catch
    display('Error: Test ut_meanMCabs_g is wrongly coded. We skip it.')
    %fprintf(fid,'Error: Test ut_meanMCabs_g is wrongly coded. We skip it.\n');
end

if 0
try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_cubMC_g);
    results=run(ut_cubMC_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
catch
    display('Error: Test ut_cubMC_g is wrongly coded. We skip it.')
    fprintf(fid,'Error: Test ut_cubMC_g is wrongly coded. We skip it.\n');
end
end

try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_cubMCabs_g);
    results=run(ut_cubMCabs_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
catch
    display('Error: Test ut_cubMCabs_g is wrongly coded. We skip it.')
    %fprintf(fid,'Error: Test ut_cubMCabs_g is wrongly coded. We skip it.\n');
end


%try
%    Tests = matlab.unittest.TestSuite.fromClass(?ut_integralNoPenalty_g);
%    results=run(ut_integralNoPenalty_g);
%    if sum([results.Failed])>0
%        failed=find([results.Failed]>0);
%        for i=1:size(failed,2)
%            fprintf(fid,'%s\n',Tests(failed(i)).Name);
%        end
%    end
%catch
%    display('Error: Test ut_integralNoPenalty_g is wrongly coded. We skip it.')
%end



time=toc;
disp(time)

diary off
