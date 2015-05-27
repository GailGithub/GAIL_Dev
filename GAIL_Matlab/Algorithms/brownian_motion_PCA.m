function Z = brownian_motion_PCA(N,T)
%% Brownian Motion by PCA
% The function brownian_motion_PCA returns $N$ $k$-steps sample of a brownian motion given the following inputs:
%
% * A vector $T=[t_1,t_2,\dots,t_k]$, with $0<t_1<t_2<\dots<t_k$ representing the points in time.
% * An integer $N$, which is the number of samples.     

Sigma=bsxfun(@min,T',T);  %Sigma(i,j) = min(T(i),T(j))

[Eigenvectors,Eigenvalues]=eig(Sigma); %Retrieving the eigenvectors and eigenvalues of the matrix Sigma

A = Eigenvectors*Eigenvalues.^(1/2);  %Defining A, which is the matrix responsible for the trasformation of the random variables into the brownian sample

Z=randn(N,length(T))*A';  %Creating N independent brownian motions

end

