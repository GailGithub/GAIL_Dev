% fasttests_cubLattice_g: fast tests for cubLattice_g

%% CALL DOCTESTS
tic; doctest cubLattice_gCLASS; time=toc


%% CALL UNIT TESTS
[~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION < 8.1 
    warning('Cannot run unit tests in MATLAB version before 8.1');
else
    try
        Tests = matlab.unittest.TestSuite.fromClass(?ut_cubLattice_gCLASS);
        results=run(ut_cubLattice_gCLASS)
        if sum([results.Failed])>0
            failed=find([results.Failed]>0);
            for i=1:size(failed,2)
                fprintf(fid,'%s\n',Tests(failed(i)).Name);
            end
        end
    catch
        display('Error: Test ut_cubLattice_gCLASS is wrongly coded. We skip it.')
        fprintf(fid,'Error: Test ut_cubLattice_g is wrongly coded. We skip it.\n');
    end
end