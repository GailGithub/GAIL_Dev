function xlat = lattice_tony(nmin,nmax, d)
% d:    dimension of the problem. Should check the maximum dimension.

y=net(sobolset(d),nmax); % Van der Corput sequence
y=y(nmin:nmax); % Which part of the sequence we need
nelem=nmax-nmin+1;

%Rank-1 lattice points
h=17797; %Korobov generator that works for a range of n=2^17
xlat=zeros(nelem,d);
xlat(:,1)=y;
for j=2:d
    xlat(:,j)=mod(xlat(:,j-1)*h,1);
end
