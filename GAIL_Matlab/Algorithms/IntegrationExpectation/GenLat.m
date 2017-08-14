function [VCool] = GenLat(FunVec)

t_start = tic;

nf=size(FunVec);
nf=nf(2); 
VNormal=[zeros(nf,2)];
V = @(x,y)

for i=1:nf
   [meanf,mean_out] = cubLattice_gCLASS(FunVec(i)); 
   for j=1:2 
      if j==1
         VNormal(i,j)=meanf-mean_out.errBd;
      elseif j==2 
         VNormal(i,j)=meanf+mean_out.errBd;
      end 
   end 
end 

VPlus= max(VNormal(:));
VMinus= min(VNormal(:));
VCool = (VPlus+VMinus)/2;


end
