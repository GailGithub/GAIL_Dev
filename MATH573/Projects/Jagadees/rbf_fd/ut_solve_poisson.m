% UT_solve_poisson unit test for solve_poisson
classdef ut_solve_poisson < matlab.unittest.TestCase
    
    methods(Test)
        % ut_solve_poisson tests solutions to the RBF-FD Approximation

        function test1(testCase)
            global GAUSSQR_PARAMETERS
            if ~isstruct(GAUSSQR_PARAMETERS)
                error('Please run test_main before running this test')
            end
            tol = 1e-4;
            g = 'squareg';
            
            [p1] = initmesh(g, 'Jiggle', 'on', 'Hmax',0.42);
            dm1 = distance_matrix(p1', p1', true);
            points = p1;
            bound_pts = union(find(points(1,:)==1), find(points(1,:)==-1));
            bound_pts = union(bound_pts, find(points(2,:)==1));
            bound_pts = union(bound_pts, find(points(2,:)==-1));
            int_pts = setdiff([1:length(points')], bound_pts);
            pi1 = int_pts; pb1 = bound_pts;
            testU = @(x,y)(sin(pi*x).*sin(pi*y));
            testLu = @(x,y)(-2.0*(pi^2)*sin(pi*x).*sin(pi*y));
            epsilon = 0.0146;
            RBFQR_flag = 1;
            [o1, o2, o3] = calc_rmse(dm1, p1,pi1,pb1, testU, testU, testLu, epsilon, RBFQR_flag);
            err_rbfdir = 0.0841; rmsed = 1.4864; stable_flag = 1;
     
            error = norm([o1, o2] - [err_rbfdir rmsed]);
            testCase.verifyLessThanOrEqual(error,tol);
            testCase.verifyEqual(o3,stable_flag);
        end

        function test2(testCase)
            global GAUSSQR_PARAMETERS
            if ~isstruct(GAUSSQR_PARAMETERS)
                error('Please run test_main before running this test')
            end
            tol = 1e-4;
            
            g = 'circleg';
            [p1] = initmesh(g, 'Jiggle', 'on', 'Hmax',0.42);
            dm1 = distance_matrix(p1', p1', true);
            points = p1;
            int_pts = find(abs(dot(points,points)-1)>1e-14);
            bound_pts = setdiff([1:length(points')], int_pts);

            pi1 = int_pts; pb1 = bound_pts;
            testU  = @(x,y) exp(-(x-0.1).^2-0.5*y.^2);
            testLu = @(x,y) exp(-(x-0.1).^2-0.5*y.^2).*(y.^2 + (-2*x+0.2).^2 - 3);
            epsilon = 0.0245;
            RBFQR_flag = 1;
            [o1, o2, o3] = calc_rmse(dm1, p1,pi1,pb1, testU, testU, testLu, epsilon, RBFQR_flag);
            
            err_rbfdir =0.0134; rmsed = 0.0942; stable_flag = 1;
            error = norm([o1, o2] - [err_rbfdir rmsed]);
            
            testCase.verifyLessThanOrEqual(error,tol);
            testCase.verifyEqual(o3,stable_flag);
       end
        
    end
end
