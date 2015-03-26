function run_handle_ut(fid, filename)
% RUN_HANDLE_UT  To run a Matlab unit tests script and handle any error
%
% Input
%   fid        Name of file .txt where results are written.
%   filename   String of Matlab filename
  try
    cls = sprintf(['?' filename]);
    %Tests = matlab.unittest.TestSuite.fromClass(cls); % This is not working
    %results=run(filename)
    eval(['Tests = matlab.unittest.TestSuite.fromClass(', cls , ')']);
    results=run(Tests)
    if sum([results.Failed])>0
      failed=find([results.Failed]>0);
      for i=1:size(failed,2)
        fprintf(fid,'%s\n',Tests(failed(i)).Name);
      end
    end
  catch
    display(['Test ' filename ' is wrongly coded. We skip it.'])
    fprintf(fid,['Test ' filename ' is wrongly coded. We skip it.\n']);
  end
end