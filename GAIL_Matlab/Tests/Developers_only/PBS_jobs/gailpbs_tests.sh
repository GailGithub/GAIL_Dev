#!/bin/sh

## Request 1 processor on 1 node
#PBS -l nodes=1:ppn=1

## Request time hh:mm:ss of walltime
#PBS -l walltime=01:00:00

## Request 1 gigabyte of memory per process
## PBS -l pmem=1gb

## set name of job
#PBS -N gail_daily_tests

## Mail alert at (b)eginning, (e)nd and (a)bortion of execution
##PBS -m bea

## Send mail to the following address
##PBS -M schoi32@iit.edu

## Use submission environment
#PBS -V

## Output files are given by cron. No need to use these.
##PBS -e gail_daily_tests.err
##PBS -o gail_daily_tests.out


####################################
# Start job from the directory it was submitted
cd $PBS_O_WORKDIR

# MATLAB
# Generate files necessary for creating HTML documentation
#./GAIL_Matlab/Documentation/Developers_only/autodoc.sh

if [ $# -eq 0 ]
  then
    echo "gailpbs_tests: No arguments supplied"
    matlabVer="R2017a"
else
  echo 'gailpbs_tests: Input arg' $1
  matlabVer=$1
fi

# Set the directory for running our matlab test
cd /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/Tests/Developers_only/automatic_tests/
# don't use this, eps files not generated  # /share/apps/matlab/R2017a/bin/matlab -nodisplay -nojvm < automatictests.m
# /share/apps/matlab/R2017a/bin/matlab -nodisplay < automatictests.m
/share/apps/matlab/$matlabVer/bin/matlab -nodisplay < automatictests.m


# SETTING THE TEST OUTPUT FILE FOR COMPARING
cp /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/OutputFiles/gail_tests-* gail_doctests.txt
mv /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/OutputFiles/gail_tests-* /home/gail/GAIL_tests/test_reports/
mv /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/OutputFiles/gail_unittests.txt gail_unittests.txt
# Preparing results txt files
begin="------Doc test wrong results were:";
nl=$'\n';
end="------Unit test wrong results were:";
echo ${begin}${nl} > test_results_grep.txt
grep -A 3 "not ok" gail_doctests.txt >> test_results_grep.txt
grep -A 3 "Warning" gail_doctests.txt | tail -n +11 >> test_results_grep.txt # We use the tail command to skip the first 2 Warnings that are regarding nojvm Matlab interface and missing admin rights to install GAIL
echo ${nl}${nl}${end}${nl} >> test_results_grep.txt
cat test_results_grep.txt gail_unittests.txt > Test_Results-$matlabVer.txt
rm test_results_grep.txt
rm gail*

# COMPARE TO THE CORRECT EXPECTED RESULTS. IF WRONG, SEND AND EMAIL WITH THE DIFFERENCES
# Compare the results obtain and send an email if different
# if diff Test_Results.txt Test_Results_correct.txt >/dev/null ; then
#  echo "Tests OK"
#  sh mail_c.sh
#else
#  echo "Tests with ERRORS"
#  sh mail_w.sh
#fi
#rm Test_Results.txt

# KEEPING ONLY THE LAST 30 DAYS REPORTS
find /home/gail/GAIL_tests/test_reports/* -mtime +30 -exec rm {} \;

# KEEPING ONLY THE LAST 30 DAYS PBS REPORTS
find /home/gail/GAIL_tests/PBS_jobs/pbs_reports/* -mtime +30 -exec rm {} \;

# send results in email
# can't call from here, output file not yet closed
# sh /home/gail/GAIL_tests/PBS_jobs/sendemail.sh $matlabVer
