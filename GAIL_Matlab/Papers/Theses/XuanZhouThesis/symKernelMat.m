function kMat = symKernelMat(kernel,dsites)
% Evaluate the kernel matrix for a symmetric kernel
N = size(dsites,1);
kMat = zeros(N);
for j = 1:N
    for k = j:N
        kMat(j,k) = kernel(dsites(j,:)',dsites(k,:)');
    end
end
kMat = kMat+triu(kMat,1)';
end