% RUN_HANDLE  To run a Matlab script and handle any runtime error
%
% Input
%   filename   String of Matlab filename
function run_handle(filename)

try
  disp('');
  disp(horzcat('Running ', filename, ' ...'));
  run(filename)
catch err
  warning(['Runtime rror executing ',filename, ':']);
  disp(err.message);
  for e=1:length(err.stack)
    stack_str = horzcat(...
      '  File: ', err.stack(e).name, ', Line ', num2str(err.stack(e).line));
    disp(stack_str);
  end
  disp('');
end