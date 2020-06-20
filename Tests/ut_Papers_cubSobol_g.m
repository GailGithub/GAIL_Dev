% UT_PAPERS_CUBSOBOL_G unit test for Papers for cubSobol_g
classdef ut_Papers_cubSobol_g < matlab.unittest.TestCase
  
  methods(Test)
    
    function cubSobol_AsianSuccessrate(testCase)
      Sobolsuccess = RunTestCubatureonGeoAsianCallSobol;
      testCase.verifyLessThanOrEqual(0.93,Sobolsuccess);
    end
    
    function cubSobol_KesiterSuccessrate(testCase)
      Sobolsuccess = RunTestCubatureonKeisterSobol;
      testCase.verifyLessThanOrEqual(0.93,Sobolsuccess);
    end
    
  end
end
