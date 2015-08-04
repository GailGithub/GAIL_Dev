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

function testKeisterDimen(testCase)
n = 3;
for d = 1:n
    for a = d/2 
        gx = @(x,b) ((sqrt(2.*pi).*b).^d).*cos(b.*sqrt(sum(x.*x,2)))...
    .*exp((1/2-b.^2).*sum(x.*x,2));
    actual = meanMC_AIS_g(gx,[0.5 2.5]);
    fQMC = @(t,a) gx(norminv(t),a);
    expected = cubSobol_g(@(t) fQMC(t,a),[zeros(1,d); ones(1,d)],'uniform', ...
      0.002,0);
  verifyEqual(testCase,actual,expected,'AbsTol',0.1)
    end
        
end


end