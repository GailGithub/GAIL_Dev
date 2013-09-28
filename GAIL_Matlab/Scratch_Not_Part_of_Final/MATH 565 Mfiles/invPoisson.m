function y=invPoisson(p,lambda)
%computes the inverse Poisson distribution
%p may be a vector, but lambda must be a scalar
if nargin<2; lambda=1; end
y=zeros(size(p));
term=exp(-lambda); %probability of exactly 0 jumps
psofar=term; %probability of at most 0 jumps
i=0; %number of jumps
notyet=1:numel(p); %which p have not been inverted yet
while numel(notyet)>0
    whgood=p(notyet)<=psofar; %these p are not inverted to i jumps
    whbad=not(whgood); %these p are too large
    y(notyet(whgood))=i; %perform inversion
    notyet=notyet(whbad); %update the p not yet inverted
    i=i+1; %update the number of jumps
    term=term*lambda/i; %update the probability of exactly i jumps
    psofar=psofar+term; %update probability of at most i jumps
end
    

