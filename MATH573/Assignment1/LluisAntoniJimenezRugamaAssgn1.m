% Lluís Antoni Jiménez Rugama, lluisantoni@gmail.com
%
% My goals for this course are 
%  - Learning the first steps needed to develop a reliable mathematical
%  software.
%  - Seeing new methods to optimize codes.
%  - Mastering Matlab.
%  - Understanding GitHub.
%  - Acknoledging better how to develop a softare in group.

tol=1e-0;
f=@(n) 1./rand(n,1); % This is not in L^1
tic
hmu=meanMC_g(f,'abstol',tol);
time=toc;
mu=Inf;
error=abs(mu-hmu);
disp(['               True solution = ' num2str(Inf)])
disp(['Our approximated MC solution = ' num2str(hmu,5)])
disp(['         Error of our method = ' num2str(error,5)])
disp(['        Tolerance guaranteed = ' num2str(tol,5)])
disp(['    Time taken for computing = ' num2str(time,5) ' seconds'])
if error>tol
    disp(['Error > tolerance!!!!!!! (O_o)'])
end

% LluisAntoniJimenezRugamaAssgn1
%                True solution = Inf
% Our approximated MC solution = 17.305
%          Error of our method = Inf
%         Tolerance guaranteed = 1
%     Time taken for computing = 0.012002 seconds
% Error > tolerance!!!!!!! (O_o)

% Problem: See that we are working with a r.v. with infinite expectaction
% and the algorithm can not detect it.

% Suggestion for improvement:
%  - One suggestion would be checking E(X^p) for a big p instead of E(X),
%  to see if it explodes really fast.
%  - Another suggestion would be checking the sample standard deviation for
%  a few n well spaced (example: comparing std(f(10^1)) and std(f(10^10)))
%  and verifying that this std_dev is increasing with probability x%.
%  instead of decreasing.
%  - The easiest solution is saying somewhere that this algorithm is only
%  working in L^1.
