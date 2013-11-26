% ut_Mutation unit test for Mutation.m
% Name: Tiago Silva
% Email: tsilva@hawk.iit.edu

classdef ut_Mutation < matlab.unittest.TestCase
    
    methods(Test)

         function UniDimensionalCaseFirtGeneration(testCase)
            tol=10^(-7);
            U=abs([0.05 0.7 0.05 ones(1,7)*0.7]);
            f = @(v,U) v+(U<0.1); T=10;
            x=0;tau_i=0; lvl=2;
            % According to the given data
            exactSol=[3 2];
            % The testSolution is
            [testSol(1,1),testSol(1,2)]=Mutation(f,x,tau_i,lvl,T,U);
            error = norm(exactSol-testSol);
            testCase.verifyLessThanOrEqual(error,tol);
         end
        
         function MutationPhaseKthGeneration(testCase)
            tol=10^(-7);
            U=abs([0.7 0.05 0.7 0.06 0.7]);
            f = @(v,U) v+(U<0.1); T=10;
            x=3;tau_i=5; lvl=5;
            % According to the given data
            exactSol=[9 5];
            % The testSolution is
            [testSol(1,1),testSol(1,2)]=Mutation(f,x,tau_i,lvl,T,U);
            error = norm(exactSol-testSol);
            testCase.verifyLessThanOrEqual(error,tol);
        end
        
        
        function MultiDimensionalCase(testCase)
            tol=10^(-7);
            U=abs(1-eye(11,10));
            f = @(v,U) v+(U<0.1); T=10;
            x=zeros(11,1);tau_i=zeros(11,1); lvl=1;
            % According to the given data
            exactSol=[(1:11)' [ones(10,1);0]];
            % The testSolution is
            testSol=zeros(11,2);
            [testSol(:,1),testSol(:,2)]=Mutation(f,x,tau_i,lvl,T,U);
            error = norm(exactSol-testSol);
            testCase.verifyLessThanOrEqual(error,tol);
        end
        
    end
    
end

% run(ut_Mutation)
% Running ut_Mutation
% ...
% Done ut_Mutation
% __________
%
%
% ans = 
%
%   1x3 TestResult array with properties:
%
%     Name
%     Passed
%     Failed
%     Incomplete
%     Duration
% 
% Totals:
%    3 Passed, 0 Failed, 0 Incomplete.
%    0.0053949 seconds testing time.