%ut_UniFunMin unit test for Xin Tong's thesis
classdef ut_UniFunMin < matlab.unittest.TestCase

  methods(Test)
      
    function testUniFunMin_test_ErrTolerance(testCase)
      nrep=10000; abstol=10^(-8); TolX=0; nmax=10^7;
      [tauvec,prob] = UniFunMin_test_ErrTolerance(nrep,abstol,TolX,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyGreaterThanOrEqual(succrates,[0.2,0.5,0.8]);
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
    function testUniFunMin_test_ErrXTolerance(testCase)
      nrep=10000; abstol=10^(-8); TolX=10^(-6); nmax=10^7;
      [tauvec,prob] = UniFunMin_test_ErrXTolerance(nrep,abstol,TolX,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyGreaterThanOrEqual(succrates,[0.2,0.5,0.8]);
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
    function testUniFunMin_test_TwoExtreme(testCase)
      nrep=10000; TolX=[10^(-2) 10^(-4) 10^(-7)]; nmax=10^7; 
      [TolXvec,prob] = UniFunMin_test_TwoExtreme(nrep,TolX,nmax);
      succrates = prob.probfunmin;   
      testCase.verifyGreaterThanOrEqual(succrates,[0.2,0.5,0.8]);
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
    function testUniFunMin_test_XTolerance(testCase)
      nrep=10000; abstol=0; TolX=10^(-6);  nmax=10^7;
      [tauvec,prob] = UniFunMin_test_XTolerance(nrep,abstol,TolX,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyGreaterThanOrEqual(succrates,[0.2,0.5,0.8]);
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
  end
end