%UT_FUNMIN_G unit test for funmin_g
classdef ut_funmin_g < matlab.unittest.TestCase

  methods(Test)
    
      
    function funmin_g1(testCase)
      f = @(x) (x-0.3).^2+1;
      in_param.abstol = 10^(-8); 
      in_param.ninit = 10;
      in_param.nmax = 10^6;
      [fmin,result] = funmin_g(f,in_param);
      realmin = 1;
      error = abs(fmin-realmin);
      testCase.verifyLessThanOrEqual(error,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funmin_g2(testCase)
      f = @(x) sin(x);
      in_param.abstol = 10^(-8); 
      in_param.ninit = 10;
      in_param.nmax = 10^6;
      [fmin,result] = funmin_g(f,in_param);
      realmin = 0;
      error = abs(fmin-realmin);
      testCase.verifyLessThanOrEqual(error,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funmin_g3(testCase)
      f = @(x) (x-pi/6).^6+1;
      in_param.abstol = 10^(-8); 
      in_param.ninit = 10;
      in_param.nmax = 10^6;
      [fmin,result] = funmin_g(f,in_param);
      realmin = 1;
      error = abs(fmin-realmin);
      testCase.verifyLessThanOrEqual(error,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
    function funmin_g4(testCase)
      f = @(x) cos(4*pi*(x-0.2));
      in_param.abstol = 10^(-8); 
      in_param.ninit = 10;
      in_param.nmax = 10^6;
      [fmin,result] = funmin_g(f,in_param);
      realmin = -1;
      error = abs(fmin-realmin);
      testCase.verifyLessThanOrEqual(error,in_param.abstol);
      testCase.verifyLessThanOrEqual(result.npoints,in_param.nmax);
    end
    
  end
end
