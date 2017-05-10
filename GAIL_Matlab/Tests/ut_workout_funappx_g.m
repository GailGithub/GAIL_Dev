%UT_WORKOUTFUNAPPX_G long unit tests for workout_funappx_g
classdef ut_workout_funappx_g < matlab.unittest.TestCase

  methods(Test)
               
%     function testworkout_timetest_funappx_g(testCase)
%       [timelgratio,~]=workout_funappx_g(100,1e-7,100,1000);
%       testCase.verifyLessThanOrEqual(timelgratio,[1.0,3.0,1.0,1.0]);
%       testCase.verifyGreaterThanOrEqual(timelgratio,[0,1.0,0,0]);
%     end
    
    function testworkout_npointstest_funappxpenalty_g(testCase)
      [~,npointslgratio]=workout_funappx_g(100,1e-7,100,1000,'funappxPenalty_g');
      testCase.verifyGreaterThanOrEqual(npointslgratio,[0.25,1.3,0.065]);
    end
    
     function testworkout_npointstest_funappx_g(testCase)
      [~,npointslgratio]=workout_funappx_g(100,1e-7,100,1000);
      testCase.verifyGreaterThanOrEqual(npointslgratio,[0.05,0.04,0.01]);
    end
  end
end
