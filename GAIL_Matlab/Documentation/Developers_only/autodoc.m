function autodoc
transform_readme;
end

function transform_readme
% find the correct directory and file
dir = which('README.txt');

% read the file
filetext = fileread(dir);

% search for title patten that consists of characters, spaces, parenthesis,
% followed by a new line of more than two hyphens. Prepend '%%' to the
% titles
title_expr = '([\w\(\) \f\v\t]+?)(\n(-{2,})\n)';
x = regexpi(filetext, title_expr, 'tokens');
replaced_text = filetext;
for i = 1:length(x),
    %[strtrim(char(x{i}(1))), '->', '% ', strtrim(char(x{i}(1))) ]
    replaced_text = regexprep(replaced_text, (char(x{i}(1))), ['\n%% ', strtrim(char(x{i}(1)))], 'once');
end
% Add '%% ' to the beginning
replaced_text = ['%% ',replaced_text];
% replace multiple hyphens in a line by itself with empty string
replaced_text = regexprep(replaced_text, '\n-{2,}', '\n');

% prepend each line with comment sign in MATLAB
replaced_text = regexprep(replaced_text, '[\n|\r|\v]', '\n% ');
% fix title lines
replaced_text = regexprep(replaced_text, '% %% ', '%% ');

% format some commands into text font in MATLAB markup langauge
algo_expr = '([a-zA-Z]+?_g(\.m){0,1})';
x3 = regexpi(replaced_text, algo_expr, 'tokens');
for i = 1:length(x),
    %[strtrim(char(x{i}(1))), '->', '% ', strtrim(char(x{i}(1))) ]
    replaced_text = regexprep(replaced_text, (char(x3{i}(1))), ['*', strtrim(char(x3{i}(1))), '*']);
end
% somehow some algorithmic names got bracketed by 2 stars---replace 
% with just 1 star
replaced_text = regexprep(replaced_text, '**', '*');
 
coms = {'doctest ', 'help '};
for i = 1:length(coms)
    replaced_text = regexprep(replaced_text, char(coms{i}), ['*',  char(coms{i}), '*']);
end

% write to help_readme.m
[GAILPATH,GAILVERSION,MATLABVERSION] = GAILstart(false);
fid = fopen([GAILPATH, 'Documentation', filesep,  'help_readme.m'],'w');
fprintf(fid, '%s', replaced_text);
fclose(fid);
end %transform_readme