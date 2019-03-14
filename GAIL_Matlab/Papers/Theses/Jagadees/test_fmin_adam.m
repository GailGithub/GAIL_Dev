function test_adam()
nDataSetSize = 1000;
vfInput = rand(1, nDataSetSize);
phiTrue = [3 2];
fhProblem = @(phi, vfInput) vfInput .* phi(1) + phi(2);
vfResp = fhProblem(phiTrue, vfInput) + randn(1, nDataSetSize) * .1;
plot(vfInput, vfResp, '.'); hold;



phi0 = randn(2, 1);
phiHat = fmin_adam(@(phi)LinearRegressionMSEGradients(phi, vfInput, vfResp), ...
  phi0, 0.01)
plot(vfInput, fhProblem(phiHat, vfInput), '.');
end


function [fMSE, vfGrad] = LinearRegressionMSEGradients(phi, vfInput, vfResp)
% - Compute mean-squared error using the current parameter estimate
vfRespHat = vfInput .* phi(1) + phi(2);
vfDiff = vfRespHat - vfResp;
fMSE = mean(vfDiff.^2) / 2;

% - Compute the gradient of MSE for each parameter
vfGrad(1) = mean(vfDiff .* vfInput);
vfGrad(2) = mean(vfDiff);
end
