function param=verifyparam(param,fieldnames,fielddim,fielddefault)
%   Make sure the parameters defining the functions exist and 
%   are the right sizes
nfield=length(fieldnames);
for i=1:nfield
    name=fieldnames{i};
    dimen=fielddim{i};
    proddim=dimen(1)*dimen(2);
    if ~isfield(param,name); param.(name)=fielddefault{i}; end
    [nx,dx]=size(param.(name));
    if (nx==dimen(1))&&(dx==dimen(2));
    elseif nx*dx >= proddim; param.(name)=reshape(param.(name)(1:proddim),dimen(1),dimen(2));
    else param.(name)=param.(name)(1)*ones(dimen(1),dimen(2));
    end
end