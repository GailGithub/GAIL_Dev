classdef DigitalOptionTest < matlab.unittest.TestCase
    % DigitalOptionTest tests the digital option price generated using
    % Monte Carlo Method and Quasi Monte Carlo Method.
    
    methods (Test)
        
        function testSobolCash(testCase)
            inp.payoffParam.optType = {'digitalcash'};
            inp.payoffParam.callPutType = {'Call'};
            inp.priceParam.cubMethod = 'Sobol';
            inp.priceParam.relTol = 0;
            inp.priceParam.absTol = 0.001;
            inp.timeDim.timeVector = 1;
            
            n = 4;
            pay = [0.8 0.9 1 1.1];
            k = [10 15 20 30];
            S = [08 18 28 38];
            i = [0.02 0.04 0.06 0.08];
            sigma = [0.4 0.45 0.5 0.55];
            
            for a = 1:n;
                for b = 1:n;
                    for c = 1:n;
                        for d = 1:n;
                            for e = 1:n;
                                inp.payoffParam.digitalPay = pay(a);
                                inp.payoffParam.strike = k(b);
                                inp.assetParam.initPrice = S(c);
                                inp.payoffParam.interest = i(d);
                                inp.assetParam.volatility = sigma(e);
                                DigitalCashCall = optPrice(inp);
                                expCashCall = DigitalCashCall.exactPrice;
                                actCashCall = genOptPrice(DigitalCashCall);
                                testCase.verifyEqual(actCashCall, ...
                                    expCashCall,'RelTol',DigitalCashCall.priceParam.relTol, ...
                                    'AbsTol',DigitalCashCall.priceParam.absTol);
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
            inp.priceParam.relTol = 0;
            inp.priceParam.absTol = 0.01;
            inp.timeDim.timeVector = 1;
            
            n = 4;
            pay = [0.8 0.9 1 1.1];
            k = [10 15 20 35];
            S = [08 18 28 38];
            i = [0.02 0.04 0.06 0.08];
            sigma = [0.4 0.45 0.5 0.55];
            
            for a = 1:n;
                for b = 1:n;
                    for c = 1:n;
                        for d = 1:n;
                            for e = 1:n;
                                inp.payoffParam.digitalPay = pay(a);
                                inp.payoffParam.strike = k(b);
                                inp.assetParam.initPrice = S(c);
                                inp.payoffParam.interest = i(d);
                                inp.assetParam.volatility = sigma(e);
                                DigitalAssetCall= optPrice(inp);
                                expAssetCall = DigitalAssetCall.exactPrice;
                                actAssetCall = genOptPrice(DigitalAssetCall);
                                testCase.verifyEqual(actAssetCall, ...
                                    expAssetCall,'RelTol',DigitalAssetCall.priceParam.relTol, ...
                                    'AbsTol',DigitalAssetCall.priceParam.absTol);
                            end
                        end
                    end
                end
            end
        end
    end
end
