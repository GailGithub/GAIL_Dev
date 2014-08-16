function nmax= estsamplebudget(tbudget,nbudget,ntry,nsofar,tstart,ttry)
% ESTSAMPLEBUDGET  estimate sample budget left from time budget by linear
% regression.
%
% Input Parameters:
%
%   tbudget---time budget
%   nbudget --- sample budget 
%   ntry --- the sample size used 
%   nsofar --- the total sample size used so far
%   tstart --- time to start the clock
%   ttry --- the time to get ntry samples
            timeleft = tbudget-toc(tstart);
            % update the time left by subtracting the time used from time
            % budget
            nleft = floor(timeleft*ntry/ttry);
            % estimate sample left by linear regression          
            nmax = max(min(nbudget-nsofar,nleft),1);
            % update the sample left 
end