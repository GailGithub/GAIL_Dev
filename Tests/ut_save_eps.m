% ut_save_eps unit tests for gail.save_eps
classdef ut_save_eps < matlab.unittest.TestCase
    
    methods(Test)
        
        function test_save_eps_1(testCase)
            GAILPATH = GAILstart(false);
            dir_name = strcat(GAILPATH,'OutputFiles',filesep,'temp',filesep);
            x = 0:0.1:1; y = sin(x);
            h = figure(1);
            plot(x,y)
            eps_file_name = gail.save_eps('temp', 'sin');
            testCase.verifyEqual(exist(eps_file_name, 'file'), 2);
            testCase.verifyEqual(length(eps_file_name), 20 + length(strcat(dir_name,'sin.eps')));
            close(h);
            rmdir(dir_name,'s');
        end
        
        function test_save_eps_2a(testCase)
            GAILPATH = GAILstart(false);
            dir_name = strcat(GAILPATH,'OutputFiles',filesep,'temp',filesep);
            x = 0:0.1:1; y = sin(x);
            h = figure(1);
            plot(x,y)
            eps_file_name = gail.save_eps('temp', 'sin', true);
            testCase.verifyEqual(exist(eps_file_name, 'file'), 2);
            testCase.verifyEqual(length(eps_file_name), 20 + length(strcat(dir_name,'sin.eps')));
            close(h);
            rmdir(dir_name,'s');
        end
        
        function test_save_eps_2b(testCase)
            GAILPATH = GAILstart(false);
            dir_name = strcat(GAILPATH,'OutputFiles',filesep,'temp',filesep);
            x = 0:0.1:1; y = sin(x);
            h = figure(1);
            plot(x,y)
            eps_file_name = gail.save_eps('temp', 'sin', false);
            testCase.verifyEqual(exist(eps_file_name, 'file'), 2);
            testCase.verifyEqual(length(eps_file_name), length(strcat(dir_name,'sin.eps')));
            close(h);
            rmdir(dir_name,'s');
        end
    end
end