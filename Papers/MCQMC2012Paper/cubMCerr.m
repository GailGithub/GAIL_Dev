function [param,Q]=cubMCerr(param,tstart)
%Handles errors in cubMC and cubMCparam
%to give an exit with information
%param.exit = 0   success
%             1   too many samples required
%             10  interval does not contain numbers
%             11  interval not 2 x d
%             12  interval is only a point in one direction
%             13  interval is infinite when measure is uniform
%             14  interval is not doubly infinite when measure is normal
%             20  importance sampling scale is not positive
if ~isfield(param,'exit'); return; end
if param.exit==0; return; end
switch param.exit
    case 1 %too many samples
        warning('GAIL:cubMC:maxreached', ['tried to evaluate at ' int2str(param.n*param.dim) ...
            ' samples,\n' ...
            'which is more than the allowed maximum of ' ...
            num2str(param.ndmax) ' samples\n']);
            return
    case 10; fprintf(2,'Error: interval must contain numbers\n');
    case 11; fprintf(2,'Error: interval must be 2 x d\n');
    case 12; fprintf(2,...
        'Error: interval must be more than a point in any coordinate direction\n');
    case 13; fprintf(2,'Error: interval must be finite when measure is uniform\n');
    case 14; fprintf(2,['Error: interval must be infinite in both directions' ...
        ' when measure is normal\n']);
    case 20; fprintf(2,'Error: importance sampling scale must be positive'\n');
end
param.Q=NaN;
Q=param.Q;
if nargin>1; param.time=toc(tstart); end

