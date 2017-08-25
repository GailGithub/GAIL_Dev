%UT_CONVTEST_FUNAPPX_G unit test for conv. test of funappx_g
classdef ut_convtest_funappx_g < matlab.unittest.TestCase

  methods(Test)
               
    function convtest_computationalcost_funappx_g(testCase)
      [npoints,errest,~,npointsglobal,errestglobal,~]=funappx_convtest;
      testCase.verifyGreaterThanOrEqual(npointsglobal,npoints);
      testCase.verifyGreaterThanOrEqual(sum(errest<=errestglobal),2);
    end
    
    function convtest_HighAccuracy_funappx_g(testCase)
       [fappx, result] = funappx_g(@(x)exp(x),'a',-2,'b',2,'ninit',15, 'abstol', 1e-12);
       x = rand(1000000,1)*(result.b-result.a)+result.a;
       actualerr = max(abs(fappx(x)-result.f(x)));
       testCase.verifyLessThanOrEqual(actualerr,result.abstol);
    end
    
  end
end

