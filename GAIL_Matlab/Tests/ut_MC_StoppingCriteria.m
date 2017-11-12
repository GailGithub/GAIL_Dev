%UT_MC_StoppingCriteria unit test for MC Stopping Criteria article
classdef ut_MC_StoppingCriteria < matlab.unittest.TestCase

  methods(Test)
             
    function testMC_StoppingCriteria_Keister(testCase)
      nrep = 1000; 
      succTable = KeisterCubatureExampleWiley(nrep);
      testCase.verifyGreaterThanOrEqual(table2array(succTable), 0.95);
    end
    
    function testMC_StoppingCriteria_Plots(testCase)
      PlotPoints;
    end
    
  end
end