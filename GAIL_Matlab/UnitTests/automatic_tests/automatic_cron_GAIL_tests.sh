#!/bin/sh

cd /home/lantoni/GAIL_tests/repo/gail-development/
# GIT
# Pulling the latest repository before testing
# /usr/local/bin/git --git-dir /home/lantoni/GAIL_tests/repo/gail-development/.git checkout .
/usr/local/bin/git --git-dir /home/lantoni/GAIL_tests/repo/gail-development/.git pull
# /usr/local/bin/git --git-dir /home/lantoni/GAIL_tests/repo/gail-development/.git checkout .

# MATLAB
# Set the directory for running our matlab test
# Run the file that installs GAIL and run the tests. The output files are in OutputFiles. We put all togehter since there is a permission not letting us install the path
cd /home/lantoni/GAIL_tests/repo/gail-development/GAIL_Matlab/UnitTests/automatic_tests/
# /export/apps/matlab/R2013a/bin/matlab -nojvm < tautomatic.m > /home/lantoni/GAIL_tests/test_results.txt Former line
/export/apps/matlab/R2013a/bin/matlab -nojvm < tautomatic.m

# SETTING THE TEST OUTPUT FILE FOR COMPARING
# Delete all the lines before first word "TAP" in our output file
# sed -i '/TAP/,$!d' test_results.txt
# Delete the lines containing "seconds testing time" from the test_results.txt file
# sed --in-place '/seconds testing time/d' test_results.txt 
# Former code above. Now below:
mv /home/lantoni/GAIL_tests/repo/gail-development/GAIL_Matlab/OutputFiles/gail_doctests.txt gail_doctests.txt
mv /home/lantoni/GAIL_tests/repo/gail-development/GAIL_Matlab/OutputFiles/gail_unittests.txt gail_unittests.txt
# Preparing results txt files
begin="------Doc test wrong results were:";
nl=$'\n';
end="------Unit test wrong results were:";
echo ${begin}${nl} > test_results_grep.txt
grep -A 3 "not ok" gail_doctests.txt >> test_results_grep.txt
echo ${nl}${nl}${end}${nl} >> test_results_grep.txt
cat test_results_grep.txt gail_unittests.txt > Test_Results.txt
rm test_results_grep.txt
rm gail*


# COMPARE TO THE CORRECT EXPECTED RESULTS. IF WRONG, SEND AND EMAIL WITH THE DIFFERENCES
# Compare the results obtain and send an email if different
if diff Test_Results.txt Test_Results_correct.txt >/dev/null ; then
  echo "Tests OK"
  sh mail_c.sh
else
  echo "Tests with ERRORS"
  sh mail_w.sh
fi
rm Test_Results.txt
