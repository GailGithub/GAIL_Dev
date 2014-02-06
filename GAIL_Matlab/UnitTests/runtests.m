% Drives all doctests and unit tests

[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
  'gail_tests-', datestr(now,'dd-mmm-yyyy-HH-MM-SS'),'.txt');
diary(filename)

tic; 
% Call doctest 
format short
doctest funappx_g
doctest funappxtau_g
doctest funappx01_g
doctest dt_funappx_g

doctest integral_g
doctest integraltau_g
doctest integralab_g
doctest dt_integral_g 

doctest meanMC_g
doctest cubMC_g

% Call unit tests
[~,~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION >= 8 
    run(ut_funappx_g)
    run(ut_funappxab_g)
    warning('off','MATLAB:integral_g:peaky')
    run(ut_integral_g)
    warning('on','MATLAB:integral_g:peaky')
    run(ut_integralab_g)
    run(ut_meanMC_g)
    run(ut_cubMC_g)
end

time=toc;
disp(time)

diary off