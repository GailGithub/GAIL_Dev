%Try scrambled Sobol

nrep=50;
barriercall=zeros(nrep,1);
for iii=1:nrep;
    OptionPrice
    barriercall(iii)=price.callN(sample.N);
end
std(barriercall)