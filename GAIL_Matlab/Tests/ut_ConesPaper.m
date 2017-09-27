%UT_CONESPAPER unit test for Cone Paper
classdef ut_ConesPaper < matlab.unittest.TestCase

  methods(Test)
             
    function testConepaper_test_integral_g(testCase)
%      nrep=10000; nmax=1e7; abstol=1e-8; %for a good table
      nrep=1000; nmax=1e7; abstol=1e-8; %for faster testing
      [succnowarn, succwarn, pfin] = conepaper_test_integral_g(nrep,nmax,abstol);
      succrates = succnowarn + succwarn;   
      testCase.verifyGreaterThanOrEqual(succrates,pfin);
      testCase.verifyGreaterThan(succrates(3), succrates(2));
      testCase.verifyGreaterThan(succrates(2), succrates(1));
    end
    
    function testConepaper_test_funappx_g(testCase)
      nrep=10000; nmax=1e7; abstol=1e-8;
      [succnowarn,succwarn] = conepaper_test_funappx_g(nrep,nmax,abstol);
      succrates = succnowarn + succwarn;   
      testCase.verifyGreaterThanOrEqual(succrates,[0.1,0.5,0.75]);
    end
  end
end