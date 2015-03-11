%ut_thesis_funmin01 unit test for funmin_g thesis 
classdef ut_thesis_funmin01_g < matlab.unittest.TestCase

  methods(Test)
      
    function test_UniFunMin_ErrTolerance(testCase)
      %nrep=10000; abstol=10^(-8);nmax=10^7; %for better tablle
      nrep=1000; abstol=10^(-8);nmax=10^7; %for faster testing
      [tauvec,prob] = UniFunMin_test_ErrTolerance(nrep,abstol,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
    function test_UniFunMin_XTolerance(testCase)
      %nrep=10000; TolX=10^(-6); nmax=10^7; %for better tablle
      nrep=1000; TolX=10^(-6); nmax=10^7; %for faster testing
      [tauvec,prob] = UniFunMin_test_XTolerance(nrep,TolX,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
    function test_UniFunMin_ErrXTolerance(testCase)
      %nrep=10000; abstol=10^(-8); TolX=10^(-6); nmax=10^7; %for better tablle
      nrep=1000; abstol=10^(-8); TolX=10^(-6); nmax=10^7; %for faster testing
      [tauvec,prob] = UniFunMin_test_ErrXTolerance(nrep,abstol,TolX,nmax);
      succrates = prob.succnowarn + prob.succwarn;   
      testCase.verifyLessThanOrEqual(succrates,[1,1,1]);
    end
    
    function test_UniFunMin_TwoExtreme(testCase)
      %nrep=10000; TolX=[10^(-2) 10^(-4) 10^(-7)]; nmax=10^7; %for better tablle 
      nrep=1000; TolX=[10^(-2) 10^(-4) 10^(-7)]; nmax=10^7; %for faster testing
      [TolXvec,prob] = UniFunMin_test_TwoExtreme(nrep,TolX,nmax);
      succrates1 = prob.probfunmin;   
      succrates2 = prob.probfminbnd;
      testCase.verifyLessThanOrEqual(succrates1,[1,1,1]);
      testCase.verifyLessThanOrEqual(succrates2,[1,1,1]);
    end
    
    function test_UniFunMin_Plot_Bump(testCase)
        a=0.03; z=0.4;
        y=UniFunMin_Plot_Bump(a,z);
    end
    
    function test_UniFunMin_Plot_TwoExtreme(testCase)
        a1=0.3; a2=0.75;
        y=UniFunMin_Plot_TwoExtreme(a1,a2);
    end
    
    function test_UniFunMin_Plot_Flat(testCase)
        a=0.5; b=1;
        y=UniFunMin_Plot_Flat(a,b);
    end
          
  end
end