function [ y ] = nufft(b,y)
%NUFFT My Fast Fourier transform with the base b and points y. The output
%is a vector
if size(y,1)<size(y,2)
    y=y';
end
mmin=log(prod(size(y)))/log(b);
    for l=0:mmin-1
       nl=b^l;
       nmminlm1=b^(mmin-l-1);
       ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
       coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)/(b*nl))';
       coefv=repmat(coef,nmminlm1,1);
       evenval=y(ptind);
       oddval=y(~ptind);
       y(ptind)=(evenval+coefv.*oddval)/b;
       y(~ptind)=(evenval-coefv.*oddval)/b;
    end
y=y';
end

