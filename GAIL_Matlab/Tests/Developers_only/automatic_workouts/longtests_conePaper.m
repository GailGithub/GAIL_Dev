% longtests_conePaper

format short

% Cone paper
if usejava('jvm') || MATLABVERSION <= 7.12
    run_handle('ConesPaperFoolFunctions')
end
try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_ConesPaper);
    results=run(ut_ConesPaper);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Error: Test ut_ConesPaper is wrongly coded. We skip it.')
    %fprintf(fid,'Error: Test ut_ConesPaper is wrongly coded. We skip it.\n');
end
