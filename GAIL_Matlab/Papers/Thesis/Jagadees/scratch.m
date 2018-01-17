% scratch pad
function scratch()
n=123;
x=(0:(1/n):1)'; 
plot(x, expsum(x, n-3))
end

function y = expsum(x, k)

h = 3;
i=1;

y = exp(2*pi*h*i*x);
for i=2:k
  y = y + exp(2*pi*h*i*x);
end
end


function test_kernel()

obj = cubMLESobol(); obj.demoKernel(1024,2,2,0.02)
obj = cubMLESobol(); obj.demoKernel(1024,2,2,0.2)
obj = cubMLESobol(); obj.demoKernel(1024,2,2,1)
obj = cubMLESobol(); obj.demoKernel(1024,2,2,2)
obj = cubMLESobol(); obj.demoKernel(1024,2,2,200)


end