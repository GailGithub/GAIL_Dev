% RUN_HANDLE  To run a Matlab script and handle any runtime error
%
% Input
%   filename   String of Matlab filename
%
% See also: run_handle_ut
%
function run_handle(filename)

try
  disp('');
  disp(horzcat('Running ', filename, ' ...'));
  run(filename)
catch err
  % warning(['Runtime error executing ',filename, ':']); % This line is
  % giving errors due to clearing the variable filename by try/catch
  disp(err.message);
  for e=1:length(err.stack)
    stack_str = horzcat(...
      '  File: ', err.stack(e).name, ', Line ', num2str(err.stack(e).line));
    disp(stack_str);
  end
  disp('');
end