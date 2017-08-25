%UT_integral_g unit test for integral_g
classdef ut_integral_g < matlab.unittest.TestCase    
    methods (Test)
       
        function testerrorfabstol(testCase)
            f=@(x) 3*x.^2;
            abstol=1e-7;
            actSolution = integral_g(f,'abstol',abstol);
            expSolution = 1;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end
        
        function testerrorf(testCase)
            abstol=1e-7;
            f=@(x) 3*x.^2;
            actSolution = integral_g(f,'abstol', abstol);
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
            actSolution = integral_g(f,'nlo',nlo,'nhi',nhi,'abstol',abstol);
            expSolution = 1;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end  
        
        function testerrorfnamx(testCase)
            f=@(x) 3*x.^2;
            abstol=1e-7;
            nmax=1e6;
            actSolution = integral_g(f,'nmax',nmax,'abstol',abstol);
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
            actSolution = integral_g(f,'abstol',abstol,'nlo',nlo,'nhi',nhi);
            expSolution = 1;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end

        
        function testerrorfabstolnmax(testCase)
            f=@(x) 3*x.^2;
            abstol=1e-7;
            nmax=1e6;
            actSolution = integral_g(f,'abstol',abstol,'nmax',nmax);
            expSolution = 1;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end
        
       
        function testerrorOfExp2x(testCase)
            f=@(x) exp(2*x);
            inparam.a=0; inparam.b=3; inparam.abstol=1e-10;
            [actSolution,out_param] = integral_g(f,inparam);
            expSolution = 201.214396746368;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            warning('off','GAIL:integral_g:exceedbudget')
            testCase.verifyThat(abs(actSolution-expSolution)*~out_param.exceedbudget,...
                IsLessThanOrEqualTo(out_param.abstol));
            testCase.verifyThat(abs(out_param.errest)*~out_param.exceedbudget,...
                IsLessThanOrEqualTo(out_param.abstol));
            warning('on','GAIL:integral_g:exceedbudget')

        end
        
        
        function testsinefunction(testCase)
            f=@(x) sin(x);
            abstol=1e-5;
            nmax=1e4;
            actSolution = integral_g(f,'nmax',nmax,'abstol',abstol);
            expSolution =  0.459697666988838;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end
        
    end
end


