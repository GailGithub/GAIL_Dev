function [pdp,disc]=projdis(design,M)
%% Beginning
[n,d]=size(design);
designvec=permute(repmat(design,[1 1 n]),[1 3 2]);
designvect=permute(designvec,[2 1 3]);
designdiff=designvec-designvect;
twopi=2*pi;
gammavec=(1:d)'/d;
Kxyjmat=zeros(n,n,d);
disc2=zeros(d,1);

%% Looping
if isfinite(M)
    lamb=1./((pi*(1:M)).^2);
    for i=1:M
        Kxyjmat=Kxyjmat+lamb(i)*cos((twopi*i)*designdiff);
    end
else
    designdiff=mod(designdiff,1);
    Kxyjmat=-designdiff.*(1-designdiff) + 1/6;
end
for k=1:d
    Kxymat=1+gammavec(k)*Kxyjmat;
    disc2(k)=-1+sum(sum(prod(Kxymat,3),2),1)/(n*n);
end

%% Cleanup
pdp2=polyfit([0; gammavec],[0; disc2],d);
pdp2=pdp2(d:-1:1);
pdp=sqrt(pdp2);
% disp(['pdp = ' num2str(pdp) ])
disc=sqrt(disc2(d));
% discb=sqrt(sum(pdp2));
% disp(['disc= ' num2str([disca discb])])

    
    
