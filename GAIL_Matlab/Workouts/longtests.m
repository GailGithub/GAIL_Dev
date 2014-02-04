% Drives all lengthy doctests and unit tests
tic; 
% Call doctest 
format short
doctest dt_meanMC_g_TrafficModel
doctest test_cubMC_g
 

% Call unit tests
[~,~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION >= 8 
    %run(ut_funappx_g)
     
end

time=toc;

disp(time)
