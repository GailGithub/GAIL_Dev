% UT_PAPERS_CUBLATTICE_G unit test for Papers for cubLattice_g
classdef ut_Papers_cubLattice_g < matlab.unittest.TestCase
  
  methods(Test)
    
    function cubLattice_AsianSuccessrate(testCase)
      Latticesuccess = RunTestCubatureonGeoAsianCallLattice;
      testCase.verifyLessThanOrEqual(0.90,Latticesuccess);
    end
    function cubLattice_KesiterSuccessrate(testCase)
      Latticesuccess = RunTestCubatureonKeisterLattice;
      testCase.verifyLessThanOrEqual(0.90,Latticesuccess);
    end
    
  end
end
