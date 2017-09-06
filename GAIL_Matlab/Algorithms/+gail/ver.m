function version = ver(toolbox_name)
% VER Return external MATLAB toolbox version
%
% >> gail.ver('GAIL')
% 2.2
% >> gail.ver('Chebfun')
% 5.4.0
%
v = ver;
n = toolbox_name;
for k = 1:length(v)
  %fprintf('%s\n', v(k).Name);
  pat = strcat('.*',toolbox_name,'.*');
  name = regexp(v(k).Name, pat, 'match', 'once');
  if ~isempty(name) && strfind(name,n),
     version = v(k).Version;
     % fprintf('%s\n', v(k).Version);
  end
end
