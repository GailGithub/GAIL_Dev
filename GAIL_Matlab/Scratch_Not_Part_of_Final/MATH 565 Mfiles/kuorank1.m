function x=kuorank1(d,n)
load kuorank1gen.mat
vander=net(sobolset(1),n);
x=mod(vander*kuogen(1:min(d,size(kuogen,1)))',1);

