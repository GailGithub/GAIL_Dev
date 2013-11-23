% DM = distance_matrix(dsites,ctrs)
% Forms the distance matrix of two sets of points in R^s,
% i.e., DM(i,j) = || datasite_i - center_j ||_2.
% Input
%   dsites: Mxs matrix representing a set of M data sites in R^s
%              (i.e., each row contains one s-dimensional point)
%   ctrs:   Nxs matrix representing a set of N centers in R^s
%              (one center per row)
% Output
%   DM:     MxN matrix whose i,j position contains the Euclidean
%              distance between the i-th data site and j-th center
  function DM = distance_matrix(dsites,ctrs,sqrt_flag)
  
  if ~exist('sqrt_flag', 'var')
      sqrt_flag = false;
  end
  M = size(dsites,1); N = size(ctrs,1);
% Algorithm is based on expanding the terms and computing each term
% explicitly, i.e.  
%         (x1 - x2)^2 = x1.^2 + x2.^2 - 2*x1*x2;
% Should be faster
  T1 = sum(dsites.*dsites,2);
  T2 = -2*dsites*ctrs';
  T3 = (sum(ctrs.*ctrs,2))';
  if sqrt_flag
	DM = sqrt(T1(:,ones(N,1)) + T2 + T3(ones(M,1),:));
  else
	DM = (T1(:,ones(N,1)) + T2 + T3(ones(M,1),:));
  end