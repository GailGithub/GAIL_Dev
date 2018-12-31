% longtests_funmin_g : long tests for funmin_g

format short

%% Workouts
% Workouts for Traub paper
run(ut_traub_paper,'test_traub_paper_funmin_g');

% funmin_g
%doctest par_funmin_g

% funmin_g
% try
%     Tests = matlab.unittest.TestSuite.fromClass(?ut_workout_funmin_g);
%     results=run(ut_workout_funmin_g);
%     if sum([results.Failed])>0
%         failed=find([results.Failed]>0);
%         %for i=1:size(failed,2)
%         %  fprintf(fid,'%s\n',Tests(failed(i)).Name);
%         %end
%     end
% catch
%     display('Error: Test ut_workout_funmin_g is wrongly coded. We skip it.')
%     %fprintf(fid,'Error: Test ut_workout_funmin_g is wrongly coded. We skip it.\n');
% end

% try
%     Tests = matlab.unittest.TestSuite.fromClass(?ut_par_funmin_g);
%     results = run(ut_par_funmin_g);
%     if sum([results.Failed])>0
%         failed=find([results.Failed]>0);
%         for i=1:size(failed,2)
%             fprintf(fid,'%s\n',Tests(failed(i)).Name);
%         end
%     end
% catch
%     display('Error: Test ut_par_funmin_g is wrongly coded. We skip it.')
% end

% try
%     Tests = matlab.unittest.TestSuite.fromClass(?ut_par_funmin_g_end);
%     results = run(ut_par_funmin_g_end);
%     if sum([results.Failed])>0
%         failed=find([results.Failed]>0);
%         for i=1:size(failed,2)
%             fprintf(fid,'%s\n',Tests(failed(i)).Name);
%         end
%     end
% catch
%     display('Error: Test ut_par_funmin_g_end is wrongly coded. We skip it.')
% end





% Function minimization thesis
% try
%     if usejava('jvm') || MATLABVERSION <= 7.12
%         run_handle('UniFunMin_Plot_Bump')
%         run_handle('UniFunMin_Plot_Flat')
%         run_handle('UniFunMin_Plot_TwoExtreme')
%     end
% catch
%     display('Error: Test for Papers/UniFunMin is wrongly coded. We skip it.')
%     %fprintf(fid,'Error: Test for Papers/UniFunMin is wrongly coded. We skip it.\n');
% end
% try
%     Tests = matlab.unittest.TestSuite.fromClass(?ut_thesis_funmin01_g);
%     results=run(ut_thesis_funmin01_g);
%     if sum([results.Failed])>0
%         failed=find([results.Failed]>0);
%         % for i=1:size(failed,2)
%         %    fprintf(fid,'%s\n',Tests(failed(i)).Name);
%         % end
%     end
% catch
%     display('Error: Test ut_thesis_funmin01_g is wrongly coded. We skip it.')
% end
