function y = fmin_ex1(x)
if exist('fmin_ex1X.mat','file')
   load fmin_ex1X xAll fAll
else
   xAll = [];
   fAll = [];
end
xAll = [xAll; x(:)];
y = -5*exp(-100*(x-0.15).^2) - exp(-80*(x-0.65).^2);
fAll = [fAll; y(:)];
save fmin_ex1X xAll fAll
end