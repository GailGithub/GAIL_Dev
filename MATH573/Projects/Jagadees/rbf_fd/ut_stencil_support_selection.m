% UT_stencil_support_selection unit test for stencil_support_selection
classdef ut_stencil_support_selection < matlab.unittest.TestCase
    
    methods(Test)
        % ut_stencil_support_selection tests solutions to the stencil_support_selection

        function test1(testCase)
			g = 'squareg';
            [p1] = initmesh(g, 'Jiggle', 'on', 'Hmax',0.42);
            dm1 = distance_matrix(p1', p1', true);	
			o1 = stencil_support_selection(dm1, p1, 25);
			stencil_support = [25     5    38    20    39     1    26];
            testCase.verifyEqual(stencil_support,o1);
        end

        function test2(testCase)
			g = 'squareg';
            [p1] = initmesh(g, 'Jiggle', 'on', 'Hmax',0.42);
            dm1 = distance_matrix(p1', p1', true);	
			o1 = stencil_support_selection(dm1, p1, 33);
			stencil_support = [33    47    43    55    26    48    39];
            testCase.verifyEqual(stencil_support,o1);
        end

        function test2(testCase)
			g = 'circleg';
			[p1,e1,t1] = initmesh(g, 'Jiggle', 'on', 'Hmax',0.42);
			[p2] = refinemesh(g,p1,e1,t1);
            dm2 = distance_matrix(p2', p2', true);
            o2 = stencil_support_selection(dm2, p2, 38);
			stencil_support = [38   130   127   129   128   126    84];
			testCase.verifyEqual(stencil_support,o2);
        end
        
    end
end

