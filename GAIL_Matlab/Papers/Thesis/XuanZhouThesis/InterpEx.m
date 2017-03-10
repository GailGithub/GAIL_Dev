% DistanceMatrixFit
% Script that uses Euclidean distance matrices to perform
% scattered data interpolation for arbitrary space dimensions
% Calls on: DistanceMatrixB, CreatePoints, testfunctionsD,
%           PlotSurf, PlotError2D, PlotSlices, PlotErrorSlices
% Uses:     various routines called by CreatePoints
kmax = 12;
rbf = @(e,r) exp(-(e*r).^2); ep = 1;
RMS_err = zeros(12,3);
N = 2.^(1:kmax);
M = 1000;
for k = 1:kmax
   for s = 1:3
      % Use Halton points as data sites and centers
      dsites = CreatePoints(N(k),s,'h');
      ctrs = dsites;
      % Create neval^s equally spaced evaluation locations in the
      % s-dimensional unit cube (should be gridded for plots)
      epoints = CreatePoints(M,s,'u');
      % Create right-hand side vector,
      % i.e., evaluate the test function at the data sites
      rhs = testfunctionsD(dsites);
      % Compute distance matrix for the data sites and centers
      DM_data = DistanceMatrixB(dsites,ctrs);
      IM = rbf(ep,DM_data);
      % Compute distance matrix for evaluation points and centers
      DM_eval = DistanceMatrixB(epoints,ctrs);
      EM = rbf(ep,DM_eval);
      % Evaluate the interpolant on evaluation points
      % (evaluation matrix * solution of interpolation system)
      Pf = EM * (IM\rhs);
      % Compute exact solution,
      % i.e., evaluate test function on evaluation points
      exact = testfunctionsD(epoints);
      % Compute maximum and RMS errors on evaluation grid
      maxerr = norm(Pf-exact,inf);
      %rms_err = norm(Pf-exact)/sqrt(M);
      RMS_err(k,s) = norm(Pf-exact)/sqrt(M);
      %fprintf('RMS error:     %e\n', rms_err)
      %fprintf('Maximum error: %e\n', maxerr)
   end
end
loglog(N,RMS_err);
xlabel('n','FontSize',12);
ylabel('RMSE','FontSize',12);
title('Interpolation Example','FontSize',14);
legend('d = 1','d = 2','d = 3');