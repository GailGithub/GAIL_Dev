% longtests_cubLattice_g : long tests for cublattice_g
format short

doctest dt_cubLattice_g

% cubLattice_g paper
try
  lattice_example;
catch
    display('Error: lattice_example is wrongly coded. We skip it.')
end
try
  FourierCoeffDecayPict;
catch
    display('Error: FourierCoeffDecayPict is wrongly coded. We skip it.')
end
%run_handle('RunTestCubatureonGeoAsianCallLattice');
%run_handle('RunTestCubatureonKeisterLattice');
try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_Papers_cubLattice_g);
    results=run(ut_Papers_cubLattice_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Error: Test ut_Papers_cubLattice_g is wrongly coded. We skip it.')
    %fprintf(fid,'Error: Test ut_Papers_cubLattice_g is wrongly coded. We skip it.\n');
end
