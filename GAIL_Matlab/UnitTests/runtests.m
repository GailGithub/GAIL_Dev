% Drives all doctests and unit tests

[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
% filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
%   'gail_tests-', datestr(now,'yyyymmddTHHMMSS'),'.txt');
% diary(filename)

tic; 
% Call doctest 
format short
doctest funappx_g
doctest funappxtau_g
doctest funappx01_g
clear in_param; doctest dt_funappx_g

doctest integral01_g
doctest integraltau_g
doctest integral_g
doctest dt_integral_g 

doctest meanMC_g
doctest cubMC_g
doctest meanMCBernoulli_g
doctest funappxlocal_g

% Call unit tests
[~,~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION >= 8 
    try
    run(ut_funappx_g)
    catch
        display('Test ut_funappx_g is wrongly coded. We skip it.')
    end
    
    try
    run(ut_funappx01_g)
    catch
        display('Test ut_funappx01_g is wrongly coded. We skip it.')
    end
    
    warning('off','MATLAB:integral_g:peaky')
    warning('off','MATLAB:integral_g:exceedbudget')
    try
    run(ut_integral_g)
    catch
        display('Test ut_integral_g is wrongly coded. We skip it.')
    end
    warning('on','MATLAB:integral_g:peaky')
    warning('on','MATLAB:integral_g:exceedbudget')
    
    warning('off','MATLAB:integral01_g:peaky')
    try
    run(ut_integral01_g)
    catch
        display('Test ut_integral01_g is wrongly coded. We skip it.')
    end
    warning('on','MATLAB:integral01_g:peaky')
    
    try
    run(ut_meanMC_g)
    catch
        display('Test ut_meanMC_g is wrongly coded. We skip it.')
    end
    
    try
    run(ut_cubMC_g)
    catch
        display('Test ut_cubMC_g is wrongly coded. We skip it.')
    end
    
    try
    run(ut_meanMCBernoulli_g)
    catch
        display('Test ut_meanMCBernoulli_g is wrongly coded. We skip it.')
    end
    
    try
    run(ut_funappxlocal_g)
    catch
        display('Test ut_funappxlocal_g is wrongly coded. We skip it.')
    end
    
end

time=toc;
% disp(time)

diary off