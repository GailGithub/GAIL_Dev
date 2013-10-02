%% Quadratic equation
clear all, close all, format long e

%% Ordinary quadratic formula
quadeq=@(a,b,c) (-b +[-1 1]'*sqrt(b^2-4*a*c))/(2*a);

%% Trying several cases
abccases=[1      3      2
          1      3e8    2
          1e200  3e200  2e200
          1e-200 3e-200 2e-200
          1      0      0
          0      1      1
          0      0      1
          0      0      0];
ncases=size(abccases,1);
    
for i=1:ncases
    a=abccases(i,1); 
    b=abccases(i,2); 
    c=abccases(i,3);
    xq=quadeq(a,b,c); %one-liner
    xr=roots([a b c]); %MATLAB's own
    xqeq=qeq(a,b,c); %an .m file
    disp('The solutions of the quadratic equation ')
    disp(['   ' num2str(a) ' x^2 + ' num2str(b) ' x + ' num2str(c) ' = 0'])
    disp('   using ''quadeq'' are:')
    fprintf('%16.14e,  %16.14e',xq)
    fprintf('\n   using ''roots'' are:\n')
    fprintf('%16.14e,  %16.14e',xr)
    fprintf('\n   using ''qeq'' are:\n')
    fprintf('%16.14e,  %16.14e',xqeq)
    fprintf('\n\n')
end
