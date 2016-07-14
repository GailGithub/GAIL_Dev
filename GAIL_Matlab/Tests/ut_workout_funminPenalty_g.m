%ut_workout_funminPenalty_g unit tests for funminPenalty_g workouts
classdef ut_workout_funminPenalty_g < matlab.unittest.TestCase

  methods(Test)
      
    function test_workout_funminPenalty_g_ErrTolerance(testCase)
      %nrep=10000; abstol=10^(-8);nmax=10^7; %for better tablle
      nrep=1000; abstol=10^(-8);nmax=10^7; %for faster testing
      [tauvec,prob] = workout_ErrToleranceTest(nrep,abstol,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
      testCase.verifyGreaterThanOrEqual(succrates,[0.3,0.4,0.7]);
    end
    
    function test_workout_funminPenalty_g_XTolerance(testCase)
      %nrep=10000; TolX=10^(-6); nmax=10^7; %for better tablle
      nrep=1000; TolX=10^(-6); nmax=10^7; %for faster testing
      [tauvec,prob] = workout_XToleranceTest(nrep,TolX,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
      testCase.verifyGreaterThanOrEqual(succrates,[0.1,0.4,0.7]);
    end
    
    function test_workout_funminPenalty_g_ErrXTolerance(testCase)
      %nrep=10000; abstol=10^(-8); TolX=10^(-6); nmax=10^7; %for better tablle
      nrep=1000; abstol=10^(-8); TolX=10^(-6); nmax=10^7; %for faster testing
      [tauvec,prob] = workout_ErrXToleranceTest(nrep,abstol,TolX,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
      testCase.verifyGreaterThanOrEqual(succrates,[0.1,0.4,0.8]);
    end
    
    function test_workout_funminPenalty_g_TwoExtreme(testCase)
      %nrep=10000; TolX=[10^(-2) 10^(-4) 10^(-7)]; nmax=10^7; %for better tablle 
      nrep=1000; TolX=[10^(-2) 10^(-4) 10^(-7)]; nmax=10^7; %for faster testing
      [TolXvec,prob] = workout_TwoExtremeTest(nrep,TolX,nmax);
      succrates1 = prob.probfunmin;   
      succrates2 = prob.probfminbnd;
      testCase.verifyLessThanOrEqual(succrates1,[1,1,1]);
      testCase.verifyLessThanOrEqual(succrates2,[1,1,1]);
      testCase.verifyGreaterThanOrEqual(succrates1,[0.6,0.6,0.6]);
      testCase.verifyGreaterThanOrEqual(succrates2,[0.6,0.6,0.6]);
    end
    
  end
end