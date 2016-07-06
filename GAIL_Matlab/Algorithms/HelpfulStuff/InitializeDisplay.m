%% InitializeDisplay: Set the display parameters to make the display beautiful
format compact %eliminate blank lines in output
close all %close all figures
set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24, ... %make font larger
      'defaultLineLineWidth',3, ... %thick lines
      'defaultLineMarkerSize',40) %big dots
LatexInterpreter %LaTeX interpreted axis labels, tick labels, and legends
%Next are MATLAB's plot colors
MATLABBlue = [0, 0.447, 0.741];
MATLABOrange = [0.85,  0.325, 0.098];
MATLABPurple = [0.494,  0.184, 0.556];
MATLABGreen = [0.466,  0.674, 0.188];
MATLABDkOrange = [0.85,  0.325, 0.098]*0.6;
MATLABLtOrange = 0.5*[0.85,  0.325, 0.098] + 0.5*[1 1 1];
