#!/bin/sh

cd /home/gail/GAIL_tests/repo/gail-development/
# GIT
# Pulling the latest repository before testing
git clean -f
/usr/local/bin/git --git-dir /home/gail/GAIL_tests/repo/gail-development/.git pull

# MATLAB
# Set the directory for running our matlab workouts
# Run the file that installs GAIL and run the tests. The output files are in OutputFiles. We put all togehter since there is a permission not letting us install the path
cd /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/Tests/Developers_only/automatic_workouts
/export/apps/matlab/R2013a/bin/matlab -nojvm < automaticworkouts.m

# SETTING THE TEST OUTPUT FILE FOR COMPARING
# Delete all the lines before first word "TAP" in our output file
# sed -i '/TAP/,$!d' test_results.txt
# Delete the lines containing "seconds testing time" from the test_results.txt file
# sed --in-place '/seconds testing time/d' test_results.txt 
# Former code above. Now below:
cd /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/OutputFiles/
find -name '*.mat' -exec mv {} /home/gail/GAIL_tests/workout_reports/ \; # Finds all the files in directories and subdirectories with extension .mat
find -name '*.eps' -exec mv {} /home/gail/GAIL_tests/workout_reports/ \; # Finds all the files in directories and subdirectories with extension .eps
# mv "/home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/OutputFiles"/*.mat /home/gail/GAIL_tests/workout_reports/
mv /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/OutputFiles/gail_workouts-* /home/gail/GAIL_tests/workout_reports/

# KEEPING ONLY THE LAST 30 DAYS REPORTS
find /home/gail/GAIL_tests/workout_reports/* -mtime +30 -exec rm {} \;
