#!/bin/sh

## Request 1 processor on 1 node
#PBS -l nodes=1:ppn=1

## Request time hh:mm:ss of walltime
#PBS -l walltime=00:01:00

## Request 1 gigabyte of memory per process
## PBS -l pmem=1gb

## set name of job
#PBS -N gail_daily_tests_email

## Mail alert at (b)eginning, (e)nd and (a)bortion of execution
##PBS -m bea

## Send mail to the following address
##PBS -M schoi32@iit.edu

## Use submission environment
#PBS -V

if [ $# -eq 0 ]
  then
    echo "Sendemail: No arguments supplied"
    matlabVer="R2016a"
else
  echo 'Sendemail: Input arg' $1
  matlabVer=$1
fi

# COMPARE TO THE CORRECT EXPECTED RESULTS. IF WRONG, SEND AND EMAIL WITH THE DIFFERENCES
# This is run after doing the daily test.
# Compare the results obtain and send an email if different
cd /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/Tests/Developers_only/automatic_tests/

if diff Test_Results-$matlabVer.txt Test_Results_correct.txt >/dev/null ; then
  echo "Tests OK"
  sh mail_c.sh $matlabVer
else
  echo "Tests with ERRORS"
  sh mail_w.sh $matlabVer
fi
rm Test_Results-$matlabVer.txt

cd -
rm gail_daily_tests_email*
