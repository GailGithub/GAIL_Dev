% ut_save_mat unit tests for gail.save_mat
classdef ut_save_mat < matlab.unittest.TestCase
    
    methods(Test)
        
        function test_save_mat_1(testCase)
            GAILPATH = GAILstart(false);
            dir_name = strcat(GAILPATH,'OutputFiles',filesep,'temp',filesep);
            x = 1; y = 2;
            mat_file_name = gail.save_mat('temp', 'x_y', true, x, y);
            clear x y;
            load(mat_file_name);
            testCase.verifyEqual(length(mat_file_name), 20 + length(strcat(dir_name,'x_y.mat')));
            testCase.verifyEqual([x,y],[1,2]);
            rmdir(dir_name,'s');
        end
        
        function test_save_mat_2(testCase)
            GAILPATH = GAILstart(false);
            dir_name = strcat(GAILPATH,'OutputFiles',filesep,'temp',filesep);
            x = 1; y = 2;
            mat_file_name = gail.save_mat('temp', 'x_y', false, x, y);
            clear x y;
            load(mat_file_name);
            testCase.verifyEqual(mat_file_name, strcat(dir_name,'x_y.mat'));
            testCase.verifyEqual([x,y],[1,2]);
            rmdir(dir_name,'s');
        end
    end
end