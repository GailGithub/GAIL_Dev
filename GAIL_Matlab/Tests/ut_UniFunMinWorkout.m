%ut_UniFunMinWorkout unit test for funmin_g Workouts
classdef ut_UniFunMinWorkout < matlab.unittest.TestCase

  methods(Test)
      
    function testWorkout_ErrToleranceTest(testCase)
      nrep=10000; abstol=10^(-8); TolX=0; nmax=10^7;
      [tauvec,prob] = workout_ErrToleranceTest(nrep,abstol,TolX,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyGreaterThanOrEqual(succrates,[0.2,0.5,0.8]);
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
    function testWorkout_ErrXToleranceTest(testCase)
      nrep=10000; abstol=10^(-8); TolX=10^(-6); nmax=10^7;
      [tauvec,prob] = workout_ErrXToleranceTest(nrep,abstol,TolX,nmax)
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyGreaterThanOrEqual(succrates,[0.2,0.5,0.8]);
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
    function testWorkout_TwoExtremeTest(testCase)
      nrep=10000; TolX=[10^(-2) 10^(-4) 10^(-7)]; nmax=10^7; 
      [TolXvec,prob] = workout_TwoExtremeTest(nrep,TolX,nmax);
      succrates = prob.probfunmin;   
      testCase.verifyGreaterThanOrEqual(succrates,[0.2,0.5,0.8]);
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
    function testWorkout_XToleranceTest(testCase)
      nrep=10000; abstol=0; TolX=10^(-6);  nmax=10^7;
      [tauvec,prob] = workout_XToleranceTest(nrep,abstol,TolX,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyGreaterThanOrEqual(succrates,[0.2,0.5,0.8]);
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
  end
end