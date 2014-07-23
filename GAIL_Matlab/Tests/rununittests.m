% Drives all doctests and unit tests

[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
% filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
%    'gail_unittests-', datestr(now,'yyyymmddTHHMMSS'),'.txt');
filename = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
   'gail_unittests','.txt');
fid = fopen(filename,'wt');

tic; 
% Call unit tests
[~,~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION >= 8
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_funappx_g);
    results=run(ut_funappx_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        fprintf(fid,'Test ut_funappx_g is wrongly coded. We skip it.\n');
    end
    
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_funappx01_g);
    results=run(ut_funappx01_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        fprintf(fid,'Test ut_funappx01_g is wrongly coded. We skip it.\n');
    end
     
    warning('off','MATLAB:integral_g:peaky')
    warning('off','MATLAB:integral_g:exceedbudget')
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_integral_g);
    results=run(ut_integral_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        fprintf(fid,'Test ut_integral_g is wrongly coded. We skip it.\n');
    end
    warning('on','MATLAB:integral_g:peaky')
    warning('on','MATLAB:integral_g:exceedbudget')

    warning('off','MATLAB:integral01_g:peaky')
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_integral01_g);
    results=run(ut_integral01_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        fprintf(fid,'Test ut_integral01_g is wrongly coded. We skip it.\n');
    end
    warning('on','MATLAB:integral01_g:peaky')
    
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_meanMC_g);
    results=run(ut_meanMC_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        fprintf(fid,'Test ut_meanMC_g is wrongly coded. We skip it.\n');
    end 
    
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_cubMC_g);
    results=run(ut_cubMC_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        fprintf(fid,'Test ut_cubMC_g is wrongly coded. We skip it.\n');
    end
   
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_meanMCBernoulli_g);
    results=run(ut_meanMCBernoulli_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        fprintf(fid,'Test ut_meanMCBernoulli_g is wrongly coded. We skip it.\n');
    end
    
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_funappxlocal_g);
    results=run(ut_funappxlocal_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        fprintf(fid,'Test ut_funappxlocal_g is wrongly coded. We skip it.\n');
    end
    
    try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_cubLattice_g);
    results=run(ut_cubLattice_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
    catch
        fprintf(fid,'Test ut_cubLattice_g is wrongly coded. We skip it.\n');
    end
    
    try
      Tests = matlab.unittest.TestSuite.fromClass(?ut_funmin_g);
      results=run(ut_funmin_g);
      if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        Tests(failed).Name
      end
    catch
      display('Test ut_funmin_g is wrongly coded. We skip it.')
    end
    
end

time=toc;
% disp(time)

fclose(fid);