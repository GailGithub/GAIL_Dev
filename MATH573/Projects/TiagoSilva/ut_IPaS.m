% ut_IPaS unit test for IPaS.m
% Name: Tiago Silva
% Email: tsilva@hawk.iit.edu

classdef ut_IPaS < matlab.unittest.TestCase
    
    methods(Test)

         function IPaSvanillaCase(testCase)
            %Naive Case using simple MC methods
            split = 3; T=10;
            f = @(v,U) v+(U<0.1);
            %std is 1% of its true value;  
            tol=0.01; 
            %true solution
            gamma = 1-binocdf(split-1,T,0.1); 
            %Sample size that satisfy solution based on MC methods
            M=ceil((1/tol)^2*(1-gamma)/gamma); 
            %Estimating gamma
            gamma_hat = IPaS(f,split,T,M);
            %Estimating variance using MC methods
            sig = sqrt(gamma_hat*(1-gamma_hat))/sqrt(M); fudge=1.2;
            error = norm(gamma-gamma_hat);
            testCase.verifyLessThanOrEqual(error,2.575*fudge*sig);
        end
        
        
        function IPaS_ExtremeCase(testCase)
            %Extreme Case
            split = [2,3,4,5,6,7]; T=10;
            f = @(v,U) v+(U<0.1); M=9000;
            % For the conditions above, the estimation of the standard 
            % deviation of the estimation is around:   
            s_hat = 1.4681e-06; 
            % So the tolerance is a 99% confidence interval
            fudge=1.2; tol= 2.575*fudge*s_hat;
            testSol = IPaS(f,split,T,M);
            exactSol = 1-binocdf(split(end)-1,10,0.1);
            error = norm(exactSol-testSol);
            testCase.verifyLessThanOrEqual(error,tol);
        end
        
           function IPaS_structured(testCase)
            %Structured Test
            coeff.split=[2,4,6]; coeff.M=9000;
            % For the conditions above, the estimation of the standard 
            % deviation of the estimation is around:  
            s_hat = 2.3279e-05; 
            % So the tolerance is a 99% confidence interval
            fudge=1.2; tol= 2.575*fudge*s_hat;
            testSol = IPaS(coeff);
            exactSol = 1-binocdf(coeff.split(end)-1,10,0.1);
            error = norm(exactSol-testSol);
            testCase.verifyLessThanOrEqual(error,tol);
        end
        
    end
    
end