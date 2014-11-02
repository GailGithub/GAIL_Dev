function newf=transformIntegrand(oldf,hyperbox,out_param)
% Transform integrand linearly so that the hyperbox would not be changed
    if strcmp(out_param.measure,'uniform') %uniform measure
        a=hyperbox(1,:); %left endpoint
        b=hyperbox(2,:); %right endpoint
    if all(a==0) && all(b==1) %no change needed
        newf=oldf; 
    else %transform points and integrand
        bmina=b-a; %hyperbox width
        volbox=prod(bmina); %volume of the hyperbox
        newf=@(x) oldf(x.*repmat(bmina,size(x,1),1)+repmat(a,size(x,1),1))...
            .*volbox;
       %stretch and shift, then multiply by volume
    end
    elseif strcmp(out_param.measure,'normal')
        newf=oldf;% no change if it is normal measure.
    end   
end