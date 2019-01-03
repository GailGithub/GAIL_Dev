function run_handle_ut(filename, varargin)
% RUN_HANDLE_UT  To run a Matlab unit tests script and handle any error
%
% Input:
%   filename   String of Matlab filename
%   testname   String of test cases to run
%   fid        Name of file .txt where results are written. Default to
%              standard output (screen)
%
% Example:
%    run_handle_ut('ut_save_eps')
%

if nargin <= 1
  fid = 1; % standard output
  testname = '';
elseif nargin <= 2
  fid = 1;
  testname = varargin{1};
else
  fid = varargin{2};
end

try
  cls = sprintf(['?' filename]);
  if isempty(testname)
    %Tests = matlab.unittest.TestSuite.fromClass(cls); % This is not working
    %results=run(filename)
    eval(['Tests = matlab.unittest.TestSuite.fromClass(', cls , ');']);
  else
    cls = sprintf(['?' filename]);
    eval(['Tests = matlab.unittest.TestSuite.fromMethod(', cls , ',''', testname, ''');']);
  end
  newline = char(13);
  disp(newline)
  results = run(Tests);
  disp('Totals:')
  disp(['   ', num2str(nnz([results.Passed])), ' Passed, ', ...
    num2str(nnz([results.Failed])), ' Failed, ', ...
    num2str(nnz([results.Incomplete])), ' Incomplete.', ...
    newline, '   ', num2str(sum([results.Duration])), ' seconds testing time.'])
  percentPassed = 100 * nnz([results.Passed]) / numel(results);
  disp(['   ', num2str(percentPassed), '% Passed.', newline]);
  rt = table(results);
  rt = sortrows(rt(:, [1:5]),'Duration');
  rt = table2cell(rt);
  [rows,~] = size(rt);
  maxlength = 40;
  for r=1:rows
    maxlength = max(length(rt{r,1}), maxlength);
  end
  minlength = maxlength;
  fprintf(1,'%*.*s %8s %8s %10s %10s\n', minlength, maxlength, 'Name', 'Passed', 'Failed', 'Incomplete', 'Duration');
  for r=1:rows
    fprintf(1,'%*.*s %8s %8s %10s     %5.4f\n', minlength, maxlength, rt{r,1}, num2str(rt{r,2}), num2str(rt{r,3}), num2str(rt{r,4}), rt{r,5});
  end
  disp(newline)
catch
  disp(['Test ', filename, ' is wrongly coded. We skip it.'])
  if nargin > 2
    fprintf(fid, '%s\n', horzcat ('Test ', filename, ' is wrongly coded. We skip it.\n'));
  end
end
end