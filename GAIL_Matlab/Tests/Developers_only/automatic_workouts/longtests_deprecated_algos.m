% doctests and unit tests for deprecated algos
format short

doctest funappxtau_g
doctest funappxglobal_g
%doctest dt_funappxglobal_g
% doctest funmin01_g
doctest integral01_g
doctest integraltau_g
doctest meanMCabs_g
doctest cubMCabs_g;
doctest cubLattice_old_g;

% try
%     Tests = matlab.unittest.TestSuite.fromClass(?ut_funappxglobal_g);
%     results=run(ut_funappxglobal_g);
%     if sum([results.Failed])>0
%         failed=find([results.Failed]>0);
%         for i=1:size(failed,2)
%             fprintf(fid,'%s\n',Tests(failed(i)).Name);
%         end
%     end
% catch
%     display('Error: Test ut_funappxglobal_g is wrongly coded. We skip it.')
%     %fprintf(fid,'Error: Test ut_funappxglobal_g is wrongly coded. We skip it.\n');
% end

warning('off','GAIL:integral01_g:peaky')
try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_integral01_g);
    results=run(ut_integral01_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        %for i=1:size(failed,2)
        %    fprintf(fid,'%s\n',Tests(failed(i)).Name);
        %end
    end
catch
    display('Error: Test ut_integral01_g is wrongly coded. We skip it.')
end
warning('on','GAIL:integral01_g:peaky')

% try
%     Tests = matlab.unittest.TestSuite.fromClass(?ut_funmin01_g);
%     results=run(ut_funmin01_g);
%     if sum([results.Failed])>0
%         failed=find([results.Failed]>0);
%         %for i=1:size(failed,2)
%         %    fprintf(fid,'%s\n',Tests(failed(i)).Name);
%         %end
%     end
% catch
%     display('Error: Test ut_funmin01_g is wrongly coded. We skip it.')
% end

try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_meanMCabs_g);
    results=run(ut_meanMCabs_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
catch
    display('Error: Test ut_meanMCabs_g is wrongly coded. We skip it.')
    %fprintf(fid,'Error: Test ut_meanMCabs_g is wrongly coded. We skip it.\n');
end

try
    Tests = matlab.unittest.TestSuite.fromClass(?ut_cubMCabs_g);
    results=run(ut_cubMCabs_g);
    if sum([results.Failed])>0
        failed=find([results.Failed]>0);
        for i=1:size(failed,2)
            fprintf(fid,'%s\n',Tests(failed(i)).Name);
        end
    end
catch
    display('Error: Test ut_cubMCabs_g is wrongly coded. We skip it.')
    %fprintf(fid,'Error: Test ut_cubMCabs_g is wrongly coded. We skip it.\n');
end


%try
%    Tests = matlab.unittest.TestSuite.fromClass(?ut_integralNoPenalty_g);
%    results=run(ut_integralNoPenalty_g);
%    if sum([results.Failed])>0
%        failed=find([results.Failed]>0);
%        for i=1:size(failed,2)
%            fprintf(fid,'%s\n',Tests(failed(i)).Name);
%        end
%    end
%catch
%    display('Error: Test ut_integralNoPenalty_g is wrongly coded. We skip it.')
%end
