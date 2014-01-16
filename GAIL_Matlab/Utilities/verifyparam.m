function in_param=verifyparam(in_param,fieldnames,fielddim,fielddefault)
%VERIFYPARAM  To make sure the parameters defining a function exist  
%   and are of the right sizes
nfield=length(fieldnames);
for i=1:nfield
  name=fieldnames{i};
  dimen=fielddim{i};
  proddim=dimen(1)*dimen(2);
  if ~isfield(in_param,name); in_param.(name)=fielddefault{i}; end
  [nx,dx]=size(in_param.(name));
  if (nx==dimen(1))&&(dx==dimen(2));
  elseif nx*dx >= proddim; in_param.(name)=reshape(in_param.(name)(1:proddim),dimen(1),dimen(2));
  else in_param.(name)=in_param.(name)(1)*ones(dimen(1),dimen(2));
  end
end