%% Brownian Motion by PCA

%%
% The function brownian_motion_PCA returns a $k$-steps sample of a brownian motion given the following inputs:
%
% * A vector $T=[t_1,t_2,\dots,t_k]$, with $0<t_1<t_2<\dots<t_k$ representing the points in time.
% * An integer $dim$, which is the number of dimensions.
function brownian = brownian_motion_PCA(T,dim)
    Sigma = zeros(length(T)); %Alocating space for the matrix Sigma, which is the matrix of covariances.
    for i=1:length(T) %Sigma(i,j) = min(T(i),T(j)) = T(min(i,j)) (The last equality is true because T(1)<T(2)<...<T(k).
        Sigma(i,i:length(T)) = T(i); 
        Sigma(i:length(T),i) = T(i);
    end
    
    [Eigenvectors,Eigenvalues]=eig(Sigma); %Retrieving the eigenvectors and eigenvalues of the matrix Sigma
    A = Eigenvectors*Eigenvalues.^(1/2);  %Defining A, which is the matrix responsible for the trasformation of the random variables into the brownian sample

    Z=A*randn(length(T),dim);  %Creating dim independent brownian motions, which each one will be a component of the brownian motion path.
    Z %Finally returns the sample
end

