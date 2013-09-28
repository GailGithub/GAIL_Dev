%Fred J. Hickernell, hickernell@iit.edu
%
%My goals for this course are 
%   o to see the students learn some important principles about writing
%     numerical software that can be used by others
%   o experience team teaching successfully, and
%   o learn some good software tricks from the others in this course.
% 

abstol=1e-12;
tic
appxinteg=integral_g(@(x) (2/sqrt(pi))*exp(-x.^2),'abstol',abstol);
time=toc;
trueinteg=erf(1);
error=trueinteg-appxinteg;
disp(['       true integral = ' num2str(trueinteg,14)])
disp(['approximate integral = ' num2str(appxinteg,14)])
disp(['               error = ' num2str(error,5)])
disp(['                 tol = ' num2str(abstol,5)])
disp(['          time taken = ' num2str(time,5) ' seconds'])

% FredHickernellAssgn1
%        true integral = 0.84270079294971
% approximate integral = 0.84270079294971
%                error = 3.2196e-15
%                  tol = 1e-12
%           time taken = 0.37054 seconds

%Suggestion for improvement: the input parsing and checkin should be moved
%below the main algorithm so that the programs are more readable.

