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

try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_Papers_cubSobol_g);
    results=run(ut_Papers_cubSobol_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
    end
catch ME
    display('Error: Test ut_Papers_cubSobol_g is wrongly coded. We skip it.')
    msgText = getReport(ME); display(msgText)
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
