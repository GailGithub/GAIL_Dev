%ut_ConesPaper unit test for Cone Paper
classdef ut_ConesPaper < matlab.unittest.TestCase

  methods(Test)
             
    function testConepaper_test_integral_g(testCase)
      conepaper_test_integral_g
      succrates = succnowarn + succwarn   
      testCase.verifyGreaterThanOrEqual(succrates,[0.1,0.5,0.8]);
    end
    
    function testConepaper_test_funappx_g(testCase)
      conepaper_test_funappx_g
      succrates = succnowarn + succwarn   
      testCase.verifyGreaterThanOrEqual(succrates,[0.1,0.5,0.8]);
    end
  end
end