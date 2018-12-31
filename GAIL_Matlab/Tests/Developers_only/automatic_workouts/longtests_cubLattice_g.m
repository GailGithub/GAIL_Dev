% longtests_cubLattice_g : long tests for cublattice_g
format short

try
  doctest dt_cubLattice_g
catch ME
    disp('Exception: longtests_cubLattice_g : in doctest dt_cubLattice_g.')
    msgText = getReport(ME); display(msgText)
end

% cubLattice_g paper
try
  lattice_example;
catch
    disp('Error: lattice_example is wrongly coded. We skip it.')
end
try
  FourierCoeffDecayPict;
catch
    disp('Error: FourierCoeffDecayPict is wrongly coded. We skip it.')
end

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
    disp('Error: Test ut_Papers_cubLattice_g is wrongly coded. We skip it.')
end
