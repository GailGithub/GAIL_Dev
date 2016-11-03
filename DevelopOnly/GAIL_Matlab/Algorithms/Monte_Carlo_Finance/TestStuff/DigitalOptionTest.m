classdef DigitalOptionTest < matlab.unittest.TestCase
    % DigitalOptionTest tests the digital option price generated using
    % Monte Carlo Method and Quasi Monte Carlo Method.
    
    methods (Test)
        
        function testSobolCash(testCase)
            inp.payoffParam.optType = {'digitalcash'};
            inp.payoffParam.callPutType = {'Call'};
            inp.priceParam.cubMethod = 'Sobol';
            inp.priceParam.relTol = 0.01;
            inp.priceParam.absTol = 0.001;
            inp.timeDim.timeVector = 1;
            
            S = 10;
            inp.assetParam.initPrice = S;
            
            n = 4;
            pay = [0.8 0.9 1 1.1];
            k = [0.7*S 0.9*S 1.1*S 1.3*S];
            i = [0.02 0.04 0.06 0.08];
            sigma = [0.4 0.45 0.5 0.55];
            
            for a = 1:n;
                for b = 1:n;
                    for c = 1:n;
                        for d = 1:n;
                            for e = 1:n;
                                inp.payoffParam.digitalPay = pay(a);
                                inp.payoffParam.strike = k(b);
                                inp.payoffParam.interest = i(c);
                                inp.assetParam.volatility = sigma(d);
                                DigitalCashCall = optPrice(inp);
                                DigitalCashPut = optPrice(inp);
                                DigitalCashPut.payoffParam.callPutType = {'Put'};
                                expCashCall = DigitalCashCall.exactPrice;
                                expCashPut = DigitalCashPut.exactPrice;
                                actCashCall = genOptPrice(DigitalCashCall);
                                actCashPut = genOptPrice(DigitalCashPut);
                                testCase.verifyEqual(actCashCall, ...
                                    expCashCall,'RelTol',DigitalCashCall.priceParam.relTol, ...
                                    'AbsTol',DigitalCashCall.priceParam.absTol);
                                testCase.verifyEqual(actCashPut, ...
                                    expCashPut,'RelTol',DigitalCashPut.priceParam.relTol, ...
                                    'AbsTol',DigitalCashPut.priceParam.absTol);
                            end
                        end
                    end
                end
            end
        end
        
        function testSobolAsset(testCase)
            inp.payoffParam.optType = {'digitalasset'};
            inp.payoffParam.callPutType = {'Call'};
            inp.priceParam.cubMethod = 'Sobol';
            inp.priceParam.relTol = 0.01;
            inp.priceParam.absTol = 0.001;
            inp.timeDim.timeVector = 1;
            
            S = 10;
            inp.assetParam.initPrice = S;
            
            n = 4;
            pay = [0.8 0.9 1 1.1];
            k = [0.7*S 0.9*S 1.1*S 1.3*S];
            i = [0.02 0.04 0.06 0.08];
            sigma = [0.4 0.45 0.5 0.55];
            
            for a = 1:n;
                for b = 1:n;
                    for c = 1:n;
                        for d = 1:n;
                                inp.payoffParam.digitalPay = pay(a);
                                inp.payoffParam.strike = k(b);
                                inp.payoffParam.interest = i(c);
                                inp.assetParam.volatility = sigma(d);
                                DigitalAssetCall = optPrice(inp);
                                DigitalAssetPut = optPrice(inp);
                                DigitalAssetPut.payoffParam.callPutType = {'Put'};
                                expAssetCall = DigitalAssetCall.exactPrice;
                                expAssetPut = DigitalAssetPut.exactPrice;
                                actAssetCall = genOptPrice(DigitalAssetCall);
                                actAssetPut = genOptPrice(DigitalAssetPut);
                                testCase.verifyEqual(actAssetCall, ...
                                    expAssetCall,'RelTol',DigitalAssetCall.priceParam.relTol, ...
                                    'AbsTol',DigitalAssetCall.priceParam.absTol);
                                testCase.verifyEqual(actAssetPut, ...
                                    expAssetPut,'RelTol',DigitalAssetPut.priceParam.relTol, ...
                                    'AbsTol',DigitalAssetPut.priceParam.absTol);
                        end
                    end
                end
            end
        end
    end
end
