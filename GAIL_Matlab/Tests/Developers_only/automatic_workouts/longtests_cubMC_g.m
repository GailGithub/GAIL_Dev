% longtests_cubMC_g: long tests for cubMC_g

% cubMC_g
format short
doctest dt_cubMC_g

% MCQMC paper
run_handle('MCQMC2012Figs')
run_handle('FoolAutomaticAlgorithms')
run_handle('RunTestcubMConGeoAsianCall')
run_handle('RunTestcubMConGaussian')
run_handle('RunTestcubMConGaussiand1')
try
    DisplayTestResults_BlacknColor({'ex1', 'ex2', 'ex3'},'black')
catch err
    display('Error: DisplayTestResults_BlacknColor is wrongly coded. We skip it.')
    display( ['file: ', err.stack(1).name, '; line number ', int2str(err.stack(1).line),  '; Error: ',  err.message ])
    %fprintf(fid,'Error: DisplayTestResults_BlacknColor is wrongly coded. We skip it.\n');
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
    display('Error: Test ut_cubMC_g is wrongly coded. We skip it.')
    fprintf(fid,'Error: Test ut_cubMC_g is wrongly coded. We skip it.\n');
end
