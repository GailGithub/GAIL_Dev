% UT_find_stencil_weights unit test for find_stencil_weights
classdef ut_find_stencil_weights < matlab.unittest.TestCase
    
    methods(Test)
        % ut_find_stencil_weights tests find_stencil_weights

        function test1(testCase)
            global GAUSSQR_PARAMETERS
            if ~isstruct(GAUSSQR_PARAMETERS)
                error('Please run test_main before running this test \n');
            end
            tol = 1e-4;
			epsilon = 0.0146;
			RBFQR_flag = 1;
			g = 'squareg';
            [p1] = initmesh(g, 'Jiggle', 'on', 'Hmax',0.42);
            dm1 = distance_matrix(p1', p1', true);	
			stencil_support = stencil_support_selection(dm1, p1, 25);
			[o1, o2, o3] = find_stencil_weights(p1(:,stencil_support)', epsilon, RBFQR_flag);
			weights = [-39.3643   10.0898    7.0783    8.1554    6.0287    7.0845    0.9276];
			v = 2.3017e-016;
			stable_flag = 1;
            error = norm(weights-o1);
            testCase.verifyLessThanOrEqual(error,tol);
			testCase.verifyLessThanOrEqual(abs(v-o2),tol);
			testCase.verifyEqual(o3,stable_flag);
        end

        function test2(testCase)
            global GAUSSQR_PARAMETERS
            if ~isstruct(GAUSSQR_PARAMETERS)
                error('Please run test_main before running this test \n');
            end
            tol = 1e-4;
			epsilon = 0.02386;
			RBFQR_flag = 1;
            
			g = 'circleg';
			[p1,e1,t1] = initmesh(g, 'Jiggle', 'on', 'Hmax',0.42);
			[p2] = refinemesh(g,p1,e1,t1);
            dm2 = distance_matrix(p2', p2', true);
			stencil_support = stencil_support_selection(dm2, p2, 38);
			[o1, o2, o3] = find_stencil_weights(p2(:,stencil_support)', epsilon, RBFQR_flag);
			weights = [-158.9434   16.9926   36.1422   34.9232   36.1027   28.3395    6.44337];
			v = 2.7395e-015;
			stable_flag = 1;
            error = norm(weights-o1);
            testCase.verifyLessThanOrEqual(error,tol);
			testCase.verifyLessThanOrEqual(abs(v-o2),tol);
			testCase.verifyEqual(o3,stable_flag);
        end
        
    end
end
