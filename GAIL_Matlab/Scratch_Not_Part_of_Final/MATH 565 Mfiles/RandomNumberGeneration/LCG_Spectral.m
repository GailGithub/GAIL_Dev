%LCG and Spectral Test
tic

M=11; %prime number
disp(['For M = ' int2str(M)])

%% Find primitive roots
aposs=(2:M-1)';
aposspow=zeros(M-2,M-2);
aposspow(:,1)=aposs;
for j=2:M-2
    aposspow(:,j)=mod(aposspow(:,j-1).*aposs,M);
end
whok=all(aposspow~=1,2);
aprim=aposs(whok);
disp('The primitive roots are')
disp(aprim')
naprim=length(aprim);

%% Perform the Spectral Test for d=2
nucand1d=(-M:M)';
[nu1cand2d,nu2cand2d]=ndgrid(nucand1d);
nucand2d=[nu1cand2d(:) nu2cand2d(:)];
nucand2d=nucand2d(any(nucand2d~=0,2),:);
dotprod=mod(nucand2d*[ones(1,naprim); aprim'],M);
whdual=dotprod==0;
l1minnorm=zeros(naprim,1);
l2minnorm=l1minnorm;
whl1minnu=zeros(naprim,2);
whl2minnu=whl1minnu;
for j=1:naprim;
    nudual2d=nucand2d(whdual(:,j),:);
    norm1dual=sum(abs(nudual2d),2);
    norm2dual=sqrt(sum(nudual2d.*nudual2d,2));
    [l1minnorm(j),whmin]=min(norm1dual);
    whl1minnu(j,:)=nudual2d(whmin(1),:);
    [l2minnorm(j),whmin]=min(norm2dual);
    whl2minnu(j,:)=nudual2d(whmin(1),:);
end
disp('The l1 spectral test is')
disp(l1minnorm')
disp('  attained at wavenumbers')
disp(whl1minnu')
disp('The l2 spectral test is')
disp(l2minnorm')
disp('  attained at wavenumbers')
disp(whl2minnu')
toc

