%% Quadratic equation
clear all, close all, format long e

%% Ordinary quadratic formula
quadeq=@(a,b,c) (-b +[-1 1]'*sqrt(b^2-4*a*c))/(2*a);

%% Roots of x^2 + 3x + 2
a=1; b=3; c=2;
xq=quadeq(a,b,c); 
xr=roots([a b c]);
xqeq=qeq(a,b,c);
disp(['The solutions of the quadratic equation ',...
    num2str(a) ' x^2 + ' num2str(b) ' x + ' num2str(c) ' = 0'])
fprintf(['   using ''quadeq'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xq)
fprintf(['    using ''roots'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xr)
fprintf(['      using ''qeq'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xqeq)
fprintf('\n')

%% Roots of x^2 + 3e8 x + 2
a=1; b=3e8; c=2;
xq=quadeq(a,b,c); 
xr=roots([a b c]);
xqeq=qeq(a,b,c);
disp(['The solutions of the quadratic equation ',...
    num2str(a) ' x^2 + ' num2str(b) ' x + ' num2str(c) ' = 0'])
fprintf(['   using ''quadeq'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xq)
fprintf(['    using ''roots'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xr)
fprintf(['      using ''qeq'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xqeq)
fprintf('\n')

%% Roots of 1e200x^2 + 3e200 x + 2e200
a=1e200; b=3e200; c=2e200;
xq=quadeq(a,b,c); 
xr=roots([a b c]);
xqeq=qeq(a,b,c);
disp(['The solutions of the quadratic equation ',...
    num2str(a) ' x^2 + ' num2str(b) ' x + ' num2str(c) ' = 0'])
fprintf(['   using ''quadeq'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xq)
fprintf(['    using ''roots'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xr)
fprintf(['      using ''qeq'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xqeq)
fprintf('\n')

%% Roots of 1e-200x^2 + 3e-200 x + 2e-200
a=1e-200; b=3e-200; c=2e-200;
xq=quadeq(a,b,c); 
xr=roots([a b c]);
xqeq=qeq(a,b,c);
disp(['The solutions of the quadratic equation ',...
    num2str(a) ' x^2 + ' num2str(b) ' x + ' num2str(c) ' = 0'])
fprintf(['   using ''quadeq'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xq)
fprintf(['    using ''roots'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xr)
fprintf(['      using ''qeq'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xqeq)
fprintf('\n')

%% Roots of x^2
a=1; b=0; c=0;
xq=quadeq(a,b,c); 
xr=roots([a b c]);
xqeq=qeq(a,b,c);
disp(['The solutions of the quadratic equation ',...
    num2str(a) ' x^2 + ' num2str(b) ' x + ' num2str(c) ' = 0'])
fprintf(['   using ''quadeq'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xq)
fprintf(['    using ''roots'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xr)
fprintf(['      using ''qeq'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xqeq)
fprintf('\n')

%% Roots of x + 1
a=0; b=1; c=1;
xq=quadeq(a,b,c); 
xr=roots([a b c]);
xqeq=qeq(a,b,c);
disp(['The solutions of the quadratic equation ',...
    num2str(a) ' x^2 + ' num2str(b) ' x + ' num2str(c) ' = 0'])
fprintf(['   using ''quadeq'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xq)
fprintf(['    using ''roots'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xr)
fprintf(['      using ''qeq'' are %16.14e, and \n' ...
    '                      %16.14e \n'],xqeq)
fprintf('\n')
