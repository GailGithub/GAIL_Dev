%ut_workoutfunappx_g unit test for workout of funappx_g
classdef ut_workoutfunappx_g < matlab.unittest.TestCase

  methods(Test)
               
    function testworkout_timetest_funappx_g(testCase)
      [timelgratio,~]=workout_funappx_g(100,1e-7,100,1000);
      testCase.verifyGreaterThanOrEqual(timelgratio,[0.4,2.0,0.1,0.75]);
    end
    
    function testworkout_npointstest_funappx_g(testCase)
      [~,npointslgratio]=workout_funappx_g(100,1e-7,100,1000);
      testCase.verifyGreaterThanOrEqual(npointslgratio,[0.25,1.3,0.065,0.3]);
    end
    
  end
end