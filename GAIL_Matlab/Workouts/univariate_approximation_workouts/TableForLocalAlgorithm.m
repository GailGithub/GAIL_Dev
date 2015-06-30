format compact
lmax=5;
npoints=zeros(lmax,4);
errorsatisfy = false(lmax,4);
nlo=10;
nhi=10;
for l=1:lmax
   l
   tol=10^(-3-l);
   [fappx,out_param,gappx,out_gparam]=testfunction(tol,nlo,nhi);
   npoints(l,1)=out_param.npoints;
   %errorsatisfy(l,1) = e  (is the 
   npoints(l,2)=out_gparam.npoints;
   [pappx1,out_param1,pappx2,out_param2]=parabolatestfunction(tol,nlo,nhi);
   npoints(l,3)=out_param1.npoints;
   npoints(l,4)=out_param2.npoints;
end
disp(npoints)
   
   