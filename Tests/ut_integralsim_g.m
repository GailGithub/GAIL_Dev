%UT_INTEGRALSIM_G unit test for integralsim_g
classdef ut_integralsim_g < matlab.unittest.TestCase    
    methods (Test)
        
        function tauchange(testCase)
            a=0.017939121772731;
            z=0.816268582603015;
            x0=z-2*a;
            x1=z+2*a;
            f = @(x) (1/(4*a^3)).*(4*a^2+(x-z).^2 ...
                +(x-z-a).*abs(x-z-a)...
                -(x-z+a).*abs(x-z+a)) ...
                .*(x>=x0).*(x<=x1); %test function
            warning('off','GAIL:integralsim_g:peaky')
            [~,out_param]=integralsim_g(f,'ninit',53,'abstol',1e-8,'nmax',1e7);
            testCase.verifyEqual(out_param.tauchange,true);
            warning('on','GAIL:integralsim_g:peaky')
        end
        
        function testerrorfabstol(testCase)
            f=@(x) 3*x.^2;
            abstol=1e-7;
            actSolution = integralsim_g(f,'abstol',abstol);
            expSolution = 1;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end
        
        function testerrorf(testCase)
            abstol=1e-7;
            f=@(x) 3*x.^2;
            actSolution = integralsim_g(f);
            expSolution = 1;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end
        
        function testerrorfninit(testCase)
            abstol=1e-7;
            f=@(x) 3*x.^2;
            ninit=53;
            actSolution = integralsim_g(f,'ninit',ninit);
            expSolution = 1;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end  
        
        function testerrorfnamx(testCase)
            f=@(x) 3*x.^2;
            abstol=1e-7;
            nmax=1e6;
            actSolution = integralsim_g(f,'nmax',nmax);
            expSolution = 1;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end

        
        function testerrorfabstolninit(testCase)
            f=@(x) 3*x.^2;
            abstol=1e-7;
            ninit=53;
            actSolution = integralsim_g(f,'abstol',abstol,'ninit',ninit);
            expSolution = 1;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end

        
        function testerrorfabstolnmax(testCase)
            f=@(x) 3*x.^2;
            abstol=1e-7;
            nmax=1e6;
            actSolution = integralsim_g(f,'abstol',abstol,'nmax',nmax);
            expSolution = 1;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end
        
        function testsinefunction(testCase)
            f=@(x) sin(x);
            abstol=1e-5;
            nmax=1e4;
            actSolution = integralsim_g(f,'nmax',nmax,'abstol',abstol);
            expSolution =  0.459697666988838;
            import matlab.unittest.constraints.IsLessThanOrEqualTo
            testCase.verifyThat(abs(actSolution-expSolution),...
                IsLessThanOrEqualTo(abstol));
        end
        
    end
end
