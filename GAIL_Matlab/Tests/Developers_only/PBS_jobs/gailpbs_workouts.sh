#!/bin/sh

## Request 1 processor on 1 node
#PBS -l nodes=1:ppn=1

## Request time hh:mm:ss of walltime
#PBS -l walltime=24:00:00

## Request 3 gigabyte of memory per process
## PBS -l pmem=1gb

## set name of job
#PBS -N gail_weekly_tests

## Mail alert at (b)eginning, (e)nd and (a)bortion of execution
##PBS -m bea

## Send mail to the following address
##PBS -M schoi32@iit.edu

## Use submission environment
#PBS -V

## Output files are given by cron. No need to use these.
##PBS -e gail_weekly_tests.err
##PBS -o gail_weekly_tests.out


####################################
# Start job from the directory it was submitted
cd $PBS_O_WORKDIR

if [ $# -eq 0 ]
  then
    echo "gailpbs_workouts: No arguments supplied"
    matlabVer="R2017a"
else
  echo 'gailpbs_workouts: Input arg' $1
  matlabVer=$1
fi

# MATLAB
# Set the directory for running our matlab workouts
# Run the file that installs GAIL and run the tests. The output files are in OutputFiles. We put all togehter since there is a permission not letting us install the path
cd /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/Tests/Developers_only/automatic_workouts/
# previous matlab version /share/apps/matlab/R2013a/bin/matlab -nojvm < automaticworkouts.m
if [ $matlabVer = "R2014a" ]
  then
    # saving figures causes segmentation fault,
    # so do not create figures
    echo 'gailpbs_workouts: disable saving figures'
    /share/apps/matlab/$matlabVer/bin/matlab -nodesktop -nojvm < automaticworkouts.m
else
  echo 'gailpbs_workouts: enable saving figures'
  /share/apps/matlab/$matlabVer/bin/matlab -nodisplay < automaticworkouts.m
fi

# SETTING THE TEST OUTPUT FILE FOR COMPARING
# Delete all the lines before first word "TAP" in our output file
# sed -i '/TAP/,$!d' test_results.txt
# Delete the lines containing "seconds testing time" from the test_results.txt file
# sed --in-place '/seconds testing time/d' test_results.txt
# Former code above. Now below:
cd /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/OutputFiles/
#find -name '*.mat' -exec mv {} /home/gail/GAIL_tests/workout_reports/ \; # Finds all the files in directories and subdirectories with extension .mat
#find -name '*.eps' -exec mv {} /home/gail/GAIL_tests/workout_reports/ \; # Finds all the files in directories and subdirectories with extension .eps

echo "Copy and then delete the eps and mat files ....."

# we use cp to retain the folder structure, whereas mv cmd does not
find -name '*.mat' -exec cp {} /home/gail/GAIL_tests/workout_reports/ \; # Finds all the files in directories and subdirectories with extension .mat
find -name '*.eps' -exec cp {} /home/gail/GAIL_tests/workout_reports/ \; # Finds all the files in directories and subdirectories with extension .eps

find -name '*.mat' -exec rm {} \; # delete it
find -name '*.eps' -exec rm {} \; # delete it

echo "Move the report files ....."

# mv "/home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/OutputFiles"/*.mat /home/gail/GAIL_tests/workout_reports/
mv /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/OutputFiles/gail_workouts-* /home/gail/GAIL_tests/workout_reports/

# Send email with the latest workout report
cd /home/gail/GAIL_tests/workout_reports/

#------------------------
WORKOUT_RESULT_FILE=$(ls -Art | grep gail_workouts- | tail -n 1)
# Preparing results txt files
begin="------Workout wrong results start:";
nl=$'\n';
end="------Workout wrong results end:";
echo ${begin}${nl} > workout_results_grep.txt
grep -A 3 "not ok" $WORKOUT_RESULT_FILE >> workout_results_grep.txt
grep -A 3 "Warning" $WORKOUT_RESULT_FILE >> workout_results_grep.txt
# | tail -n +11 We use the tail command to skip the first 2 Warnings that are regarding nojvm Matlab interface and missing admin rights to install GAIL
echo ${nl}${nl}${end}${nl} >> workout_results_grep.txt
cat workout_results_grep.txt > Workout_Results-$matlabVer.txt
rm workout_results_grep.txt
#------------------------

echo "Finally, send the report in email ....."

#mail -s $(ls -Art | grep gail_workouts- | tail -n 1) schoi32@iit.edu,jrathin1@iit.edu < $(ls -Art | grep gail_workouts- | tail -n 1)
# sending as attachment
TO_EMAIL_IDS="Sou-cheng<schoi32@iit.edu>, Jagadeeswaran<jrathin1@iit.edu>"
MAIL_BODY="Please find the workouts report attached"
SUBJECT="$matlabVer: $(ls -Art | grep gail_workouts- | tail -n 1)"

# dont send email from here
echo $MAIL_BODY | mail -s "${SUBJECT}" -a $(ls -Art | grep gail_workouts- | tail -n 1) -a Workout_Results-$matlabVer.txt ${TO_EMAIL_IDS}

# KEEPING ONLY THE LAST 30 DAYS REPORTS
find /home/gail/GAIL_tests/workout_reports/* -mtime +30 -exec rm {} \;

echo "Finished workouts !!!"
