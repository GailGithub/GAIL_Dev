% ut_save_mat unit tests for gail.save_mat 
classdef ut_save_mat < matlab.unittest.TestCase

  methods(Test)
      
    function test_save_mat_1(testCase)
      x = 1; y = 2; 
      mat_file_name = gail.save_mat('temp', 'x_y', x, y);
      clear x y;
      load(mat_file_name);
      testCase.verifyEqual([x,y],[1,2]);
      delete(mat_file_name);
      rmdir('temp')
    end
    
  end
end