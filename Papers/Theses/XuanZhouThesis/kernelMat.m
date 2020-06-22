function kMat = kernelMat(kernel,dsites,ctrs)
% Evaluate the kernel matrix for a symmetric kernel
M = size(dsites,1);
N = size(ctrs,1);
kMat = zeros(M,N);
for j = 1:M
    for k = 1:N
        kMat(j,k) = kernel(dsites(j,:)',ctrs(k,:)');
    end
end
end