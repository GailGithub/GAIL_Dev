% Drives all doctests and unit tests
format short

GAILPATH = GAILstart(0);
completereport = strcat(GAILPATH,'OutputFiles',filesep,...
  'gail_tests-', datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.txt');
diary(completereport)

shortutestreport = strcat(GAILPATH,'OutputFiles',filesep,...
   'gail_unittests','.txt');
fid = fopen(shortutestreport,'wt');

tic 
%% CALL DOCTESTS 
tic; doctest gail.gail_in_param; time=toc
tic; doctest gail.gail1D_in_param; time=toc
tic; doctest gail.funappx_g_in_param; time=toc
tic; doctest gail.funmin_g_in_param; time=toc
tic; doctest gail.integral_g_in_param; time=toc
tic; doctest gail.gailMD_in_param; time=toc
tic; doctest gail.cubMC_g_in_param; time=toc

tic; doctest funappx_g; time=toc
tic; doctest funappxNoPenalty_g; time=toc
tic; doctest dt_funappx_g; time=toc
tic; doctest dt_funappxNoPenalty_g; time=toc
tic; doctest funmin_g; time=toc
tic; doctest integral_g; time=toc
tic; doctest integralsim_g; time=toc
tic; doctest dt_integral_g ; time=toc
tic; doctest meanMC_g; time=toc
tic; doctest cubLattice_g; time=toc
tic; doctest cubSobol_g; time=toc
tic; doctest dt_integral_g; time=toc

%% CALL UNIT TESTS
[~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION < 8.1
    warning('Cannot run unit tests in MATLAB version before 8.1');
else
	warning ('off','all'); % Disable warnings. Many algorithms have warnings that will be too tedious to print
	MyFolderInfo = dir(strcat(GAILPATH,'Tests'))
	for i = 1:size(MyFolderInfo,1)
		if size(MyFolderInfo(i).name,2) > 3 % We check whether we have at least 4 characters in the name string
			if all([all(MyFolderInfo(i).name(1:3) == 'ut_') all(MyFolderInfo(i).name(4) ~= 'w')]) % We do not want ut_workouts
				run_handle_ut(MyFolderInfo(i).name(1:end-2),fid) % We substract the extension .m in name
			end
		end
	end
	warning ('on','all');
end

time=toc
% disp(time)

diary off
fclose(fid);
format
