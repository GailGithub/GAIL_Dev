%ut_workout_integral_g unit test for workout_integral_g
classdef ut_workout_integral_g < matlab.unittest.TestCase

  methods(Test)
             
    function testWorkout_integral_g(testCase)
      nrep=100; nmax=1e7; abstol=1e-8;
      [succnowarn, succwarn, pfin] = workout_integral_g(nrep,nmax,abstol);
      succrates = succnowarn + succwarn;   
      testCase.verifyGreaterThanOrEqual(succrates,pfin);
      for i=2:len(succrates)-1
        testCase.verifyGreaterThan(succrates(i), succrates(i-1));
      end
    end
    
  end
end