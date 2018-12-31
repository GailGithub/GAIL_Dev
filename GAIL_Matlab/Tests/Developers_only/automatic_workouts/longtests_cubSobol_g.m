%% longtests_cubSobol_g: long tests for cubSobol_g

format short
doctest dt_cubSobol_g

% cubSobol_g paper
try
  SobolWalshPict;
catch ME
    display('Error: SobolWalshPict is wrongly coded. We skip it.')
    msgText = getReport(ME); display(msgText)
end

try
  WalshFourierCoeffDecay;
catch ME
    display('Error: WalshFourierCoeffDecay is wrongly coded. We skip it.')
    msgText = getReport(ME); display(msgText)
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
catch ME
    display('Error: Test ut_Papers_cubSobol_g is wrongly coded. We skip it.')
    msgText = getReport(ME); display(msgText)
    %fprintf(fid,'Error: Test ut_Papers_cubSobol_g is wrongly coded. We skip it.\n');
end

% Wiley paper
try
  Tests = matlab.unittest.TestSuite.fromClass(?ut_MC_StoppingCriteria);
    results=run(ut_MC_StoppingCriteria);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
    end
catch
    display('Error: Test ut_MC_StoppingCriteria is wrongly coded. We skip it.')
end
