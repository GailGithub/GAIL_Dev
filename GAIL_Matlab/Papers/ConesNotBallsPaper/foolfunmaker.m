function y=foolfunmaker(x,piecefun,coef,c)
ncoef=length(coef);
y=zeros(size(x));
for i=1:ncoef
    y=y+coef(i)*piecefun(x,c(i,:));
end
end
