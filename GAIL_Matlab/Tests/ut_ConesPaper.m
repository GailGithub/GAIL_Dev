%ut_ConesPaper unit test for Cone Paper
classdef ut_ConesPaper < matlab.unittest.TestCase

  methods(Test)
             
    function testConepaper_test_integral_g(testCase)
      conepaper_test_integral_g
      succrates = succnowarn + succwarn   
      testCase.verifyGreaterThanOrEqual(succrates,[0.1,0.6,0.8]);
    end
    
  end
end