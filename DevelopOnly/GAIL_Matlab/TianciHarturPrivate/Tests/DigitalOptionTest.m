classdef DigitalOptionTest < matlab.unittest.TestCase
    % BasketOptionTest tests the basket option price generated using
    % Monte Carlo Method and Quasi Monte Carlo Method.

    methods (Test)
        
        function testIIDCash(testCase)
            inp.payoffParam.optType = {'digitalcash'};
            inp.payoffParam.digitalPay = 1;
            inp.assetParam.initPrice = 10;
            inp.payoffParam.interest = 0.01;
            inp.assetParam.volatility = 0.5;
            inp.priceParam.absTol = 0;
            inp.priceParam.relTol = 0.005;
            inp.timeDim.timeVector = 1;
            DigitalCashCall = optPrice(inp);
            expCashCall = DigitalCashCall.exactPrice;
            actCashCall = genOptPrice(DigitalCashCall);
            DigitalCashCall.payoffParam.putCallType = {'put'};
            DigitalCashPut = DigitalCashCall;
            expCashPut = DigitalCashPut.exactPrice;
            actCashPut = genOptPrice(DigitalCashPut);
            testCase.verifyEqual(actCashCall,expCashCall,'RelTol',DigitalCashCall.priceParam.relTol);
            testCase.verifyEqual(actCashPut,expCashPut,'RelTol',DigitalCashPut.priceParam.relTol);
        end
        
        function testSobolCash(testCase)
            inp.payoffParam.optType = {'digitalcash'};
            inp.payoffParam.digitalPay = 1;
            inp.assetParam.initPrice = 10;
            inp.assetParam.volatility = 0.5;
            inp.payoffParam.interest = 0.01;
            inp.priceParam.absTol = 0;
            inp.priceParam.relTol = 0.005;
            inp.timeDim.timeVector = 1;
            inp.priceParam.cubMethod = 'Sobol';
            DigitalCashCall = optPrice(inp);
            expCashCall = DigitalCashCall.exactPrice;
            actCashCall = genOptPrice(DigitalCashCall);
            DigitalCashCall.payoffParam.putCallType = {'put'};
            DigitalCashPut = DigitalCashCall;
            expCashPut = DigitalCashPut.exactPrice;
            actCashPut = genOptPrice(DigitalCashPut);
            testCase.verifyEqual(actCashCall,expCashCall,'RelTol',DigitalCashCall.priceParam.relTol);
            testCase.verifyEqual(actCashPut,expCashPut,'RelTol',DigitalCashPut.priceParam.relTol);
        end
        
        function testIIDAsset(testCase)
            inp.payoffParam.optType = {'digitalasset'};
            inp.payoffParam.digitalPay = 1;
            inp.assetParam.initPrice = 10;
            inp.assetParam.volatility = 0.5;
            inp.payoffParam.interest = 0.01;
            inp.priceParam.absTol = 0;
            inp.priceParam.relTol = 0.005;
            inp.timeDim.timeVector = 1;
            DigitalAssetCall = optPrice(inp);
            expAssetCall = DigitalAssetCall.exactPrice;
            actAssetCall = genOptPrice(DigitalAssetCall);
            DigitalAssetCall.payoffParam.putCallType = {'put'};
            DigitalAssetPut = DigitalAssetCall;
            expAssetPut = DigitalAssetPut.exactPrice;
            actAssetPut = genOptPrice(DigitalAssetPut);
            testCase.verifyEqual(actAssetCall,expAssetCall,'RelTol',DigitalAssetCall.priceParam.relTol);
            testCase.verifyEqual(actAssetPut,expAssetPut,'RelTol',DigitalAssetPut.priceParam.relTol);
        end
        
        function testSobolAsset(testCase)
            inp.payoffParam.optType = {'digitalasset'};
            inp.payoffParam.digitalPay = 1;
            inp.assetParam.initPrice = 10;
            inp.assetParam.volatility = 0.5;
            inp.payoffParam.interest = 0.01;
            inp.priceParam.absTol = 0;
            inp.priceParam.relTol = 0.005;
            inp.timeDim.timeVector = 1;
            inp.priceParam.cubMethod = 'Sobol';
            DigitalAssetCall = optPrice(inp);
            expAssetCall = DigitalAssetCall.exactPrice;
            actAssetCall = genOptPrice(DigitalAssetCall);
            DigitalAssetCall.payoffParam.putCallType = {'put'};
            DigitalAssetPut = DigitalAssetCall;
            expAssetPut = DigitalAssetPut.exactPrice;
            actAssetPut = genOptPrice(DigitalAssetPut);
            testCase.verifyEqual(actAssetCall,expAssetCall,'RelTol',DigitalAssetCall.priceParam.relTol);
            testCase.verifyEqual(actAssetPut,expAssetPut,'RelTol',DigitalAssetPut.priceParam.relTol);
        end
    end
end


