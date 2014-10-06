#!/bin/sh

cd /home/lantoni/GAIL_tests/repo/gail-development/
# GIT
# Pulling the latest repository before testing
/usr/local/bin/git --git-dir /home/lantoni/GAIL_tests/repo/gail-development/.git pull

# MATLAB
# Set the directory for running our matlab workouts
# Run the file that installs GAIL and run the tests. The output files are in OutputFiles. We put all togehter since there is a permission not letting us install the path
cd /home/lantoni/GAIL_tests/repo/gail-development/GAIL_Matlab/Workouts/automatic_workouts
/export/apps/matlab/R2013a/bin/matlab -nojvm < automaticworkouts.m

# SETTING THE TEST OUTPUT FILE FOR COMPARING
# Delete all the lines before first word "TAP" in our output file
# sed -i '/TAP/,$!d' test_results.txt
# Delete the lines containing "seconds testing time" from the test_results.txt file
# sed --in-place '/seconds testing time/d' test_results.txt 
# Former code above. Now below:
mv /home/lantoni/GAIL_tests/repo/gail-development/GAIL_Matlab/OutputFiles/*.mat /home/lantoni/GAIL_tests/workout_reports/

# KEEPING ONLY THE LAST 4*10 REPORTS
find /home/lantoni/GAIL_tests/workout_reports/*.mat -mtime +40 -exec rm {} \;
