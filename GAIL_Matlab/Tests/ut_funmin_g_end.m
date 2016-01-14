%UT_FUNMIN_G_end unit test for funmin_g at end points
classdef ut_funmin_g_end < matlab.unittest.TestCase
    
    methods(Test)
           
        function funmin_gEXM1(testCase)
            % minimum appears at left-end point
            f=@(x) sin(x);
            [fmin,out]=funmin_g(f);
            actualminima=0;
            testCase.verifyEqual(actualminima,out.intervals(1,1));
        end
        
        function funmin_gEXM2(testCase)
            % minimum appears at left-end point
            f=@(x) sin(1-x);
            [fmin,out]=funmin_g(f);
            actualminima=1;
            testCase.verifyEqual(actualminima,out.intervals(2,end));
        end
     
    end
end

