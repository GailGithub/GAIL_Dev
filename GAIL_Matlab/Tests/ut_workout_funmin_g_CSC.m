%ut_workout_funmin_g_CSC unit tests for funmin_g_CSC workouts
classdef ut_workout_funmin_g_CSC < matlab.unittest.TestCase

  methods(Test)
      
    function test_workout_funmin_g_CSC_ErrTolerance(testCase)
      %nrep=10000; abstol=10^(-8);nmax=10^7; %for better tablle
      nrep=1000; abstol=10^(-8);nmax=10^7; %for faster testing
      [tauvec,prob] = workout_ErrToleranceTest_CSC(nrep,abstol,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
    function test_workout_funmin_g_CSC_XTolerance(testCase)
      %nrep=10000; TolX=10^(-6); nmax=10^7; %for better tablle
      nrep=1000; TolX=10^(-6); nmax=10^7; %for faster testing
      [tauvec,prob] = workout_XToleranceTest_CSC(nrep,TolX,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
    function test_workout_funmin_g_CSC_ErrXTolerance(testCase)
      %nrep=1000; abstol=10^(-8); TolX=10^(-6); nmax=10^7; %for better tablle
      nrep=1000; abstol=10^(-8); TolX=10^(-6); nmax=10^7; %for faster testing
      [tauvec,prob] = workout_ErrXToleranceTest_CSC(nrep,abstol,TolX,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
    function test_workout_funmin_g_CSC_TwoExtreme(testCase)
      %nrep=10000; TolX=[10^(-2) 10^(-4) 10^(-7)]; nmax=10^7; %for better tablle 
      nrep=1000; TolX=[10^(-2) 10^(-4) 10^(-7)]; nmax=10^7; %for faster testing
      [TolXvec,prob] = workout_TwoExtremeTest_CSC(nrep,TolX,nmax);
      succrates1 = prob.probfunmin;   
      succrates2 = prob.probfminbnd;
      testCase.verifyLessThanOrEqual(succrates1,[1,1,1]);
      testCase.verifyLessThanOrEqual(succrates2,[1,1,1]);
    end
    
  end
end