% Drives all doctests and unit tests
 
% Call doctest 
format short
doctest funappx_g
doctest dt_funappx_g
doctest integral_g
doctest meanMC_g
% doctest cubMC_g

% Call unit tests
[~,~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION >= 8 
    run(ut_funappx_g)
    warning('off','MATLAB:integral_g:peaky')
    run(ut_integral_g)
    warning('on','MATLAB:integral_g:peaky')
    run(ut_meanMC_g)
%     run(ut_cubMC_g)
end


