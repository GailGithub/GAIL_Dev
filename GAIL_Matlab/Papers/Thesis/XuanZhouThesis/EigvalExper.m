% Experiments with eigenvalues
format compact, format short e
n=20;
tauvec=2.^(0:0.2:2)'; nt=size(tauvec,1);

v=randn(n,n);
v=v./repmat(sqrt(sum(v.*v)),n,1);

l=sort(2*abs(randn(n,1)),1,'descend');
cuml=cumsum(l);

a=v*diag(l)*v';

eigval=sort(eig(a),1,'descend');
cumeig=cumsum(eigval);

disp(' ')
disp('Eigenvalues and Diagonal Values')
disp([eigval l])

disp(' ')
disp('Cumulative Sums of Eigenvalues and Diagonal Values')
disp([cumeig cuml cumeig-cuml])

Mtau=zeros(nt,2);
for i=1:nt
    tau=tauvec(i);
    Mtau(i,:)=sum([eigval l].^(1/tau),1).^tau;
end
disp(' ')
disp('Sums of Powers of Eigenvalues and Diagonal Values')
disp([tauvec Mtau Mtau(:,1)-Mtau(:,2)])