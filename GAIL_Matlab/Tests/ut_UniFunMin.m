%ut_UniFunMin unit test for Xin Tong's thesis
classdef ut_UniFunMin < matlab.unittest.TestCase

  methods(Test)
      
    function testUniFunMin_test_ErrorTolerance(testCase)
      UniFunMin_test_ErrorTolerance
      succrates = succnowarn + succwarn;   
      testCase.verifyGreaterThanOrEqual(succrates,[0.2,0.5,0.8]);
    end
    
    function testUniFunMin_test_ErrorXTolerance(testCase)
      UniFunMin_test_ErrorXTolerance
      succrates = succnowarn + succwarn;   
      testCase.verifyGreaterThanOrEqual(succrates,[0.2,0.5,0.8]);
    end
    
    function testUniFunMin_test_TwoExtreme(testCase)
      UniFunMin_test_TwoExtreme
      succrates = probfunmin;   
      testCase.verifyGreaterThanOrEqual(succrates,[0.2,0.5,0.8]);
    end
    
    function testUniFunMin_test_XTolerance(testCase)
      UniFunMin_test_XTolerance
      succrates = succnowarn + succwarn;   
      testCase.verifyGreaterThanOrEqual(succrates,[0.2,0.5,0.8]);
    end
    
  end
end