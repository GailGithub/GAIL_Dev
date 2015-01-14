% Drives all doctests and unit tests
format short

[GAILPATH,~,PATHNAMESEPARATOR] = GAILstart(0);
completereport = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
  'gail_tests-', datestr(now,'yyyymmdd'),'.txt');
diary(completereport)

shortutestreport = strcat(GAILPATH,'OutputFiles',PATHNAMESEPARATOR,...
   'gail_unittests','.txt');
fid = fopen(shortutestreport,'wt');

tic; 
% CALL DOCTESTS 

doctest funappxglobal_g
doctest funappx_g
doctest funappxtau_g
doctest funappx01_g
doctest dt_funappxglobal_g
doctest dt_funappx_g
doctest funmin_g
doctest funmin01_g
doctest integral_g
doctest integralsim_g
doctest integral01_g
doctest integraltau_g
doctest dt_integral_g 
doctest meanMCabs_g
doctest meanMC_g
doctest meanMCBer_g
doctest cubMCabs_g
doctest cubLattice_g
doctest cubSobol_g
doctest cubLattice_old_g
doctest cubSobol_old_g

% CALL UNIT TESTS
[~,~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION >= 8
  try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_funappxglobal_g);
    results=run(ut_funappxglobal_g)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display('Test ut_funappxglobal_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_funappxglobal_g is wrongly coded. We skip it.\n');
  end
  
  try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_funappx01_g);
    results=run(ut_funappx01_g)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display('Test ut_funappx01_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_funappx01_g is wrongly coded. We skip it.\n');
  end
  
  warning('off','MATLAB:integral_g:peaky')
  warning('off','MATLAB:integral_g:exceedbudget')
  try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_integral_g);
    results=run(ut_integral_g)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display('Test ut_integral_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_integral_g is wrongly coded. We skip it.\n');
  end
  warning('on','MATLAB:integral_g:peaky')
  warning('on','MATLAB:integral_g:exceedbudget')
  
  warning('off','MATLAB:integralsim_g:peaky')
  warning('off','MATLAB:integralsim_g:exceedbudget')
  try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_integralsim_g);
    results=run(ut_integralsim_g)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display('Test ut_integralsim_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_integralsim_g is wrongly coded. We skip it.\n');
  end
  warning('on','MATLAB:integralsim_g:peaky')
  warning('on','MATLAB:integralsim_g:exceedbudget')
  
  
  warning('off','MATLAB:integral01_g:peaky')
  try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_integral01_g);
    results=run(ut_integral01_g)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display('Test ut_integral01_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_integral01_g is wrongly coded. We skip it.\n');
  end
  warning('on','MATLAB:integral01_g:peaky')
  
  try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_meanMCabs_g);
    results=run(ut_meanMCabs_g)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display('Test ut_meanMCabs_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_meanMCabs_g is wrongly coded. We skip it.\n');
  end
  
  try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_cubMCabs_g);
    results=run(ut_cubMCabs_g)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display('Test ut_cubMCabs_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_cubMCabs_g is wrongly coded. We skip it.\n');
  end
  
  try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_meanMCBer_g);
    results=run(ut_meanMCBer_g)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display('Test ut_meanMCBer_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_meanMCBer_g is wrongly coded. We skip it.\n');
  end
  
  try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_funappx_g);
    results=run(ut_funappx_g)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display('Test ut_funappx_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_funappx_g is wrongly coded. We skip it.\n');
  end
  
  try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_cubLattice_g);
    results=run(ut_cubLattice_g)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display('Test ut_cubLattice_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_cubLattice_g is wrongly coded. We skip it.\n');
  end
  
  try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_cubSobol_g);
    results=run(ut_cubSobol_g)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display('Test ut_cubSobol_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_cubSobol_g is wrongly coded. We skip it.\n');
  end
  
  
  try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_funmin_g);
    results=run(ut_funmin_g)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display('Test ut_funmin_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_funmin_g is wrongly coded. We skip it.\n');
  end
  
  try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_funmin01_g);
    results=run(ut_funmin01_g)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display('Test ut_funmin01_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_funmin01_g is wrongly coded. We skip it.\n');
  end 
  
  try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_meanMC_g);
    results=run(ut_meanMC_g)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display('Test ut_meanMC_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_meanMC_g is wrongly coded. We skip it.\n');
  end
  
  try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_cubMC_g);
    results=run(ut_cubMC_g)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display('Test utcubMC_g is wrongly coded. We skip it.')
    fprintf(fid,'Test ut_cubMC_g is wrongly coded. We skip it.\n');
  end
end

time=toc;
% disp(time)

diary off
fclose(fid);
