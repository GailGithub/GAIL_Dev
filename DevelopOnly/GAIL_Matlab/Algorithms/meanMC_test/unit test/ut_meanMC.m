%% ut_meanMC
% fast unit tests for the meanMC class
% Author: Tianpei Qian
classdef ut_meanMC < matlab.unittest.TestCase

  methods(Test)
             
    function meanMCNoInput(testCase)
       test = meanMC;
       mu = genMu(test);
       testCase.verifyClass(test,'meanMC');
       testCase.verifyLessThanOrEqual(abs(mu-1/3),0.1*1/3);
    end
       
    function plainMethod(testCase)
        obj.in_param.abstol = 2e-2;
        obj.in_param.reltol = 2e-1;
        obj.in_param.alpha = 2e-2;
        obj.in_param.fudge = 1.21;
        obj.in_param.nSig = 1.1e4;
        obj.in_param.n1 = 1.1e4;
        obj.in_param.tbudget = 200;
        obj.in_param.nbudget = 2e9;
        obj.method = {'plain'};
        obj.Yrand = @(n) sin(rand(n, 1));
        test = meanMC(obj);
        mu = genMu(test);
        testCase.verifyClass(test,'meanMC');
        testCase.verifyEqual(test.in_param.abstol, 2e-2);
        testCase.verifyEqual(test.in_param.reltol, 2e-1);
        testCase.verifyEqual(test.in_param.alpha, 2e-2);
        testCase.verifyEqual(test.in_param.fudge, 1.21);
        testCase.verifyEqual(test.in_param.nSig, 1.1e4);
        testCase.verifyEqual(test.in_param.n1, 1.1e4);
        testCase.verifyEqual(test.in_param.tbudget, 200);
        testCase.verifyEqual(test.in_param.nbudget, 2e9);
        testCase.verifyEqual(test.method, {'plain'});
        testCase.verifyLessThanOrEqual(abs(mu-(1-cos(1))),2e-2);
    end

    function cvMethod(testCase)
        obj.method = {'cv'};
        obj.cv_param.YXrand = @sinr_cv;
        obj.cv_param.muX = 0.5;
        obj.cv_param.ncv = 1.1e3;
        obj.cv_param.ridge = 1e-3;
        test = meanMC(obj);
        mu = genMu(test);
        testCase.verifyClass(test,'meanMC');
        testCase.verifyEqual(test.method, {'cv'});
        testCase.verifyEqual(test.cv_param.muX, 0.5);
        testCase.verifyEqual(test.cv_param.ncv, 1.1e3);
        testCase.verifyEqual(test.cv_param.ridge, 1e-3);
        testCase.verifyLessThanOrEqual(abs(mu-(1-cos(1))),1e-2);
    end
    
    function avMethod(testCase)
        obj.method = {'av'};
        obj.av_param.YYrand = @sinr_av;
        test = meanMC(obj);
        mu = genMu(test);
        testCase.verifyClass(test,'meanMC');
        testCase.verifyEqual(test.method, {'av'});
        testCase.verifyLessThanOrEqual(abs(mu-(1-cos(1))),1e-2);
    end
    
    function combinedMethod(testCase)
        obj.method = {'cv','av','plain'};
        obj.nc = 1e3;
        obj.Yrand = @(n) sin(rand(n, 1));
        obj.cv_param.YXrand = @sinr_cv;
        obj.cv_param.muX = 0.5;
        obj.cv_param.ridge = [1e-3 2e-3];
        obj.av_param.YYrand = @sinr_av;
        test = meanMC(obj);
        mu = genMu(test);
        testCase.verifyClass(test,'meanMC');
        testCase.verifyEqual(test.nc, 1e3);
        testCase.verifyEqual(test.method, {'cv','av','plain'});
        testCase.verifyEqual(test.cv_param.ridge, [1e-3 2e-3]);
        testCase.verifyLessThanOrEqual(abs(mu-(1-cos(1))),1e-2);
    end
    
    function multipleMus1(testCase)
        obj.in_param.abstol = [1e-2 2e-2];
        obj.in_param.reltol = [0 2e-1];
        obj.in_param.alpha = [1e-2 2e-2];
        obj.in_param.fudge = [1.21 1.22];
        obj.in_param.nSig = [1.1e4 1.2e4];
        obj.in_param.n1 = [1.1e4 1.2e4];
        obj.in_param.tbudget = [110 200];
        obj.in_param.nbudget = [1.1e9 2e9];
        obj.method = {{'cv','av','plain'},{'cv','av'}};
        obj.nc = [1e3 2e3];
        obj.Yrand = {@(n) sin(rand(n, 1)), @(n) exp(rand(n, 1))};
        obj.cv_param.YXrand = {@sinr_cv, @expr_cv};
        obj.cv_param.muX = {0.5, 0.5};
        obj.cv_param.ncv = [1e3 2e3];
        obj.cv_param.ridge = {0, [0 1e-3]};
        obj.av_param.YYrand = {@sinr_av, @expr_av};
        test = meanMC(obj);
        mu = genMu(test);
        testCase.verifyClass(test,'meanMC');
        testCase.verifyEqual(test(1).in_param.abstol, 1e-2);
        testCase.verifyEqual(test(2).in_param.abstol, 2e-2);
        testCase.verifyEqual(test(1).in_param.reltol, 0);
        testCase.verifyEqual(test(2).in_param.reltol, 2e-1);
        testCase.verifyEqual(test(1).in_param.alpha, 1e-2);
        testCase.verifyEqual(test(2).in_param.alpha, 2e-2);
        testCase.verifyEqual(test(1).in_param.fudge, 1.21);
        testCase.verifyEqual(test(2).in_param.fudge, 1.22);
        testCase.verifyEqual(test(1).in_param.nSig, 1.1e4);
        testCase.verifyEqual(test(2).in_param.nSig, 1.2e4);
        testCase.verifyEqual(test(1).in_param.n1, 1.1e4);
        testCase.verifyEqual(test(2).in_param.n1, 1.2e4);
        testCase.verifyEqual(test(1).in_param.tbudget, 110);
        testCase.verifyEqual(test(2).in_param.tbudget, 200);
        testCase.verifyEqual(test(1).in_param.nbudget, 1.1e9);
        testCase.verifyEqual(test(2).in_param.nbudget, 2e9);
        testCase.verifyEqual(test(1).method, {'cv','av','plain'});
        testCase.verifyEqual(test(2).method, {'cv','av'});
        testCase.verifyEqual(test(1).nc, 1e3);
        testCase.verifyEqual(test(2).nc, 2e3);
        testCase.verifyEqual(test(1).cv_param.muX, 0.5);
        testCase.verifyEqual(test(2).cv_param.muX, 0.5);
        testCase.verifyEqual(test(1).cv_param.ncv, 1e3);
        testCase.verifyEqual(test(2).cv_param.ncv, 2e3);
        testCase.verifyEqual(test(1).cv_param.ridge, 0);
        testCase.verifyEqual(test(2).cv_param.ridge, [0 1e-3]);
        testCase.verifyLessThanOrEqual([abs(mu(1)-(1-cos(1))), ...
            abs(mu(2)-(exp(1)-1))],[1e-2 2e-2]);
    end
    
    function multipleMus2(testCase)
        obj.in_param.abstol = {1e-2 2e-2};
        obj.in_param.reltol = {0 2e-1};
        obj.in_param.alpha = {1e-2 2e-2};
        obj.in_param.fudge = {1.21 1.22};
        obj.in_param.nSig = {1.1e4 1.2e4};
        obj.in_param.n1 = {1.1e4 1.2e4};
        obj.in_param.tbudget = {110 200};
        obj.in_param.nbudget = {1.1e9 2e9};
        obj.method = {{'cv','av','plain'},{'cv','av'}};
        obj.nc = {1e3 2e3};
        obj.Yrand = {@(n) sin(rand(n, 1)), @(n) exp(rand(n, 1))};
        obj.cv_param.YXrand = {@sinr_cv, @expr_cv};
        obj.cv_param.muX = {0.5, 0.5};
        obj.cv_param.ncv = {1e3 2e3};
        obj.cv_param.ridge = {0, [0 1e-3]};
        obj.av_param.YYrand = {@sinr_av, @expr_av};
        test = meanMC(obj);
        mu = genMu(test);
        testCase.verifyClass(test,'meanMC');
        testCase.verifyEqual(test(1).in_param.abstol, 1e-2);
        testCase.verifyEqual(test(2).in_param.abstol, 2e-2);
        testCase.verifyEqual(test(1).in_param.reltol, 0);
        testCase.verifyEqual(test(2).in_param.reltol, 2e-1);
        testCase.verifyEqual(test(1).in_param.alpha, 1e-2);
        testCase.verifyEqual(test(2).in_param.alpha, 2e-2);
        testCase.verifyEqual(test(1).in_param.fudge, 1.21);
        testCase.verifyEqual(test(2).in_param.fudge, 1.22);
        testCase.verifyEqual(test(1).in_param.nSig, 1.1e4);
        testCase.verifyEqual(test(2).in_param.nSig, 1.2e4);
        testCase.verifyEqual(test(1).in_param.n1, 1.1e4);
        testCase.verifyEqual(test(2).in_param.n1, 1.2e4);
        testCase.verifyEqual(test(1).in_param.tbudget, 110);
        testCase.verifyEqual(test(2).in_param.tbudget, 200);
        testCase.verifyEqual(test(1).in_param.nbudget, 1.1e9);
        testCase.verifyEqual(test(2).in_param.nbudget, 2e9);
        testCase.verifyEqual(test(1).method, {'cv','av','plain'});
        testCase.verifyEqual(test(2).method, {'cv','av'});
        testCase.verifyEqual(test(1).nc, 1e3);
        testCase.verifyEqual(test(2).nc, 2e3);
        testCase.verifyEqual(test(1).cv_param.muX, 0.5);
        testCase.verifyEqual(test(2).cv_param.muX, 0.5);
        testCase.verifyEqual(test(1).cv_param.ncv, 1e3);
        testCase.verifyEqual(test(2).cv_param.ncv, 2e3);
        testCase.verifyEqual(test(1).cv_param.ridge, 0);
        testCase.verifyEqual(test(2).cv_param.ridge, [0 1e-3]);
        testCase.verifyLessThanOrEqual([abs(mu(1)-(1-cos(1))), ...
            abs(mu(2)-(exp(1)-1))],[1e-2 2e-2]);
    end
    
    function multipleMus3(testCase)
        obj.in_param.abstol = 1e-2;
        obj.in_param.reltol = 0;
        obj.in_param.alpha = 1e-2;
        obj.in_param.fudge = 1.21;
        obj.in_param.nSig = 1.1e4;
        obj.in_param.n1 = 1.1e4;
        obj.in_param.tbudget = 110;
        obj.in_param.nbudget = 1.1e9;
        obj.method = {'cv','av','plain'};
        obj.nc = 1e3;
        obj.Yrand = {@(n) sin(rand(n, 1)), @(n) exp(rand(n, 1))};
        obj.cv_param.YXrand = {@sinr_cv, @expr_cv};
        obj.cv_param.muX = 0.5;
        obj.cv_param.ncv = 1e3;
        obj.cv_param.ridge = [0 1e-3];
        obj.av_param.YYrand = {@sinr_av, @expr_av};
        test = meanMC(obj);
        mu = genMu(test);
        testCase.verifyClass(test,'meanMC');
        testCase.verifyEqual(test(1).in_param.abstol, 1e-2);
        testCase.verifyEqual(test(2).in_param.abstol, 1e-2);
        testCase.verifyEqual(test(1).in_param.reltol, 0);
        testCase.verifyEqual(test(2).in_param.reltol, 0);
        testCase.verifyEqual(test(1).in_param.alpha, 1e-2);
        testCase.verifyEqual(test(2).in_param.alpha, 1e-2);
        testCase.verifyEqual(test(1).in_param.fudge, 1.21);
        testCase.verifyEqual(test(2).in_param.fudge, 1.21);
        testCase.verifyEqual(test(1).in_param.nSig, 1.1e4);
        testCase.verifyEqual(test(2).in_param.nSig, 1.1e4);
        testCase.verifyEqual(test(1).in_param.n1, 1.1e4);
        testCase.verifyEqual(test(2).in_param.n1, 1.1e4);
        testCase.verifyEqual(test(1).in_param.tbudget, 110);
        testCase.verifyEqual(test(2).in_param.tbudget, 110);
        testCase.verifyEqual(test(1).in_param.nbudget, 1.1e9);
        testCase.verifyEqual(test(2).in_param.nbudget, 1.1e9);
        testCase.verifyEqual(test(1).method, {'cv','av','plain'});
        testCase.verifyEqual(test(2).method, {'cv','av','plain'});
        testCase.verifyEqual(test(1).nc, 1e3);
        testCase.verifyEqual(test(2).nc, 1e3);
        testCase.verifyEqual(test(1).cv_param.muX, 0.5);
        testCase.verifyEqual(test(2).cv_param.muX, 0.5);
        testCase.verifyEqual(test(1).cv_param.ncv, 1e3);
        testCase.verifyEqual(test(2).cv_param.ncv, 1e3);
        testCase.verifyEqual(test(1).cv_param.ridge, [0 1e-3]);
        testCase.verifyEqual(test(2).cv_param.ridge, [0 1e-3]);
        testCase.verifyLessThanOrEqual([abs(mu(1)-(1-cos(1))), ...
            abs(mu(2)-(exp(1)-1))],1e-2);
    end
   
  end
end
