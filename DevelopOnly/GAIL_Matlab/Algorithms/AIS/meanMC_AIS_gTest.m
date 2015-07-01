function tests = meanMC_AIS_gTest
tests = functiontests(localfunctions);
end

function testMaximumExample(testCase)
actual = meanMC_AIS_g(@(t,b_value)max(t+b_value,0).*exp(-t.*b_value-...
    b_value.^2/2),[-2 2],1,0.002,0.01,1e4,1.2);
expected = 1/sqrt(2*pi);
verifyEqual(testCase,actual,expected,'AbsTol',0.01)
end

function testKeister(testCase)
actual = meanMC_AIS_g(@(x,b) ((sqrt(2.*pi).*b).^3).*cos(b.*sqrt(sum(x.*x,2)))...
    .*exp((1/2-b.^2).*sum(x.*x,2)),[0.5 3],1,0.002,0.01,1e4,1.2);
expected = 2.1686;
verifyEqual(testCase,actual,expected,'AbsTol',0.01)
end

%function testOptions(testCase)
%K = [10 15 20];
%So = [15 20 25];
%r = 0;
%d = [1 2 3];
%T = [1 2 3];
%deltaT = T./d;
%n = 3;

%SK=@(z,bval) max(0,K-mean(((So./d).*exp(cumsum((r-((sigma.^2)./2)).*deltaT + ...
%    sigma.*sqrt(deltaT).*(z-bval),2))),2));
%fx=@(z,bval) SK(z,bval).*exp(-r.*T).*exp(sum(z,2).*bval-(d.*(bval.^2)./2));

%for K(1):K(n)
%    for So(1):So(n)
%        for d(1):d(n)
%            for T(1):T(n)
%                actual = meanMC_AIS_g(fx,[0.5 3],1,0.002,0.01,1e4,1.2);
%                expected = 2.1686;
%                verifyEqual(testCase,actual,expected,'AbsTol',0.01)
%            end
%        end
%    end
%end 
%end