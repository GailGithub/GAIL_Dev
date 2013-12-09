%UT_INTEGRAL_G unit test for integral_g
classdef ut_integralab_g < matlab.unittest.TestCase    
    methods (Test)
        

        
        function testerrorfabstol(testCase)
            f=@(x) 3*x.^2;
            abstol=1e-7;
            actSolution = integralab_g(f,'abstol',abstol);
            expSolution = 1;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end
        
        function testerrorf(testCase)
            abstol=1e-7;
            f=@(x) 3*x.^2;
            actSolution = integralab_g(f);
            expSolution = 1;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end
        
        function testerrorfninit(testCase)
            abstol=1e-7;
            f=@(x) 3*x.^2;
            nlo=52;
            nhi=52;
            actSolution = integralab_g(f,'nlo',nlo,'nhi',nhi);
            expSolution = 1;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end  
        
        function testerrorfnamx(testCase)
            f=@(x) 3*x.^2;
            abstol=1e-7;
            nmax=1e6;
            actSolution = integralab_g(f,'nmax',nmax);
            expSolution = 1;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end

        
        function testerrorfabstolninit(testCase)
            f=@(x) 3*x.^2;
            abstol=1e-7;
            nlo=52;
            nhi=52;
            actSolution = integralab_g(f,'abstol',abstol,'nlo',nlo,'nhi',nhi);
            expSolution = 1;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end

        
        function testerrorfabstolnmax(testCase)
            f=@(x) 3*x.^2;
            abstol=1e-7;
            nmax=1e6;
            actSolution = integralab_g(f,'abstol',abstol,'nmax',nmax);
            expSolution = 1;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end
        
    end
end
