%ut_workout_funmin_g unit tests for funmin_g workouts
classdef ut_workout_funmin_g < matlab.unittest.TestCase

  methods(Test)
      
    function test_workout_funmin_g_ErrTolerance(testCase)
      %nrep=10000; abstol=10^(-8); nmax=10^7; %for better table
      nrep=1000; abstol=10^(-8); nmax=10^7; %for faster testing
      [tauvec,prob] = workout_ErrToleranceTest(nrep,abstol,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
    function test_workout_funmin_g_XTolerance(testCase)
      %nrep=10000; TolX=10^(-6); nmax=10^7; %for better table
      nrep=1000; TolX=10^(-6); nmax=10^7; %for faster testing
      [tauvec,prob] = workout_XToleranceTest(nrep,TolX,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
    function test_workout_funmin_g_ErrXTolerance(testCase)
      %nrep=10000; abstol=10^(-8); TolX=10^(-6); nmax=10^7; %for better table
      nrep=1000; abstol=10^(-8); TolX=10^(-6); nmax=10^7; %for faster testing
      [tauvec,prob] = workout_ErrXToleranceTest(nrep,abstol,TolX,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
    function test_workout_funmin_g_TwoExtreme(testCase)
      %nrep=10000; TolX=[10^(-2) 10^(-4) 10^(-7)]; nmax=10^7; %for better table
      nrep=1000; TolX=[10^(-2) 10^(-4) 10^(-7)]; nmax=10^7; %for faster testing
      [TolXvec,prob] = workout_TwoExtremeTest(nrep,TolX,nmax);
      succrates1 = prob.probfunmin;   
      succrates2 = prob.probfminbnd;
      testCase.verifyLessThanOrEqual(succrates1,[1,1,1]);
      testCase.verifyLessThanOrEqual(succrates2,[1,1,1]);
    end
    
  end
end